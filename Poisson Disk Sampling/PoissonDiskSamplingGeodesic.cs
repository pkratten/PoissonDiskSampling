using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;
using System.Windows.Forms;
using Grasshopper.Kernel;
using Grasshopper.Kernel.Geometry.Delaunay;
using Rhino.Display;
using Rhino.Geometry;

namespace Poisson_Disk_Sampling
{
    public class PoissonDiskSamplingGeodesic : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the PoissonDiskSamplingGeodesic class.
        /// </summary>
        public PoissonDiskSamplingGeodesic()
          : base("Poisson Disk Sampling Geodesic", "Poisson Geodesic",
              "Sample points on a surface with the poisson disk sampling method.",
              "Vector", "Grid")
        {
        }

        public override GH_Exposure Exposure => GH_Exposure.tertiary | GH_Exposure.obscure;

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            Param.Brep = pManager.AddBrepParameter("Brep", "B", "The object to sample points on.", GH_ParamAccess.item);
            Param.Distance = pManager.AddNumberParameter("Distance", "D", "The distasnce between points or the diameter of a disk.", GH_ParamAccess.item);
            Param.K = pManager.AddIntegerParameter("Number of tries", "K", "Number of tries to sample a point per cell", GH_ParamAccess.item, 3);
            Param.I = pManager.AddNumberParameter("Initial sample modifier", "I", "Modifier for the density of initial uniform random sample of points.", GH_ParamAccess.item, 1);
            Param.P = pManager.AddIntegerParameter("Number of phase groups per dimension", "P", "Number of phase groups per dimension. So the total phase groups will be P³. The minimum is 3.", GH_ParamAccess.item, 3);
            Param.Seed = pManager.AddIntegerParameter("Seed", "S", "The seed that for the random generation of the points.", GH_ParamAccess.item, 0);

            pManager[2].Optional = true;
            pManager[3].Optional = true;
            pManager[4].Optional = true;
            pManager[5].Optional = true;
        }

        struct Param
        {
            public static int Brep, Distance, K, I, P, Seed , Samples, Normals, FaceID, U, V;
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            Param.Samples = pManager.AddPointParameter("Sampled Points", "P", "Points sampled by the poisson disk sampling method.", GH_ParamAccess.list);
            Param.Normals = pManager.AddVectorParameter("Normals", "N", "Normals of the sampled points.", GH_ParamAccess.list);
            Param.FaceID = pManager.AddIntegerParameter("Face indeces", "I", "Brep faces indeces of the sampled points.", GH_ParamAccess.list);
            Param.U = pManager.AddNumberParameter("U parameter", "U", "U parameter of sampled points.", GH_ParamAccess.list);
            Param.V = pManager.AddNumberParameter("V parameter", "V", "V parameter of sampled points.", GH_ParamAccess.list);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            Brep brep = null;
            double distance = 0;
            int K = 0;
            double I = 0;
            int P = 0;
            int seed = 0;

            if (!DA.GetData(Param.Brep, ref brep)) return;
            if (!DA.GetData(Param.Distance, ref distance)) return;
            if (!DA.GetData(Param.K, ref K)) return;
            if (!DA.GetData(Param.I, ref I)) return;
            if (!DA.GetData(Param.P, ref P)) return;
            if (!DA.GetData(Param.Seed, ref seed)) return;

            if(distance <= 0)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Distance need to be positive!");
                return;
            }
            if(K < 1)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "K needs to be at least 1!");
                return;
            }
            if(I <= 0)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "I needs to be larger than 0!");
                return;
            }
            if(P < 3)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "P need to be at least 3!");
                return;
            }

            List<Point> result = new List<Point>();
            try
            {
                result = PoissonDiskSampleGeodesic(brep, distance, K, I, P, seed);
            }
            catch(Exception e)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, e.Message);
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, e.StackTrace);
                return;
            }

            List<Point3d> samples = new List<Point3d>();
            List<Vector3d> normals = new List<Vector3d>();
            List<int> faceIDs = new List<int>();
            List<double> us = new List<double>();
            List<double> vs = new List<double>();

            foreach (var point in result)
            {
                samples.Add(point.Position);
                normals.Add(point.normal);
                faceIDs.Add(point.FaceID);
                us.Add(point.PositionUV.X);
                vs.Add(point.PositionUV.Y);
            }

            DA.SetDataList(Param.Samples, samples);
            DA.SetDataList(Param.Normals, normals);
            DA.SetDataList(Param.FaceID, faceIDs);
            DA.SetDataList(Param.U, us);
            DA.SetDataList(Param.V, vs);
        }

        public class Point
        {
            public Point3d Position;
            public int FaceID;
            public Point2d PositionUV;
            public Vector3d normal;
            public Cell Cell = new Cell();
        }

        public class Cell
        {
            public int X, Y, Z;
            public Box Box;
            public int PhaseGroup = -1;
            public List<Point> InitialPoints = new List<Point>();
            public List<Point> Samples = new List<Point>();
        }

        public List<Point> PoissonDiskSampleGeodesic(Brep brep, double distance, int K, double I, int P, int seed)
        {
            Random random = new Random(seed);



            //Sample Si initial points
            ConcurrentBag<Point> initialPoints = new ConcurrentBag<Point>();
            List<Task> tasks = new List<Task>();
            for (int i = 0; i < brep.Faces.Count; i++)
            {
                tasks.Add(new Task(new Action<object>((iObj) =>
                {
                    int index = (int)iObj;

                    BrepFace face = brep.Faces[index];
                    double area = AreaMassProperties.Compute(brep.Faces[index]).Area;
                    int n = Convert.ToInt32((area / distance) * I);
                    for (int j = 0; j < n; j++)
                    {
                        Point point = new Point();
                        point.PositionUV = new Point2d(face.Domain(0).ParameterAt(random.NextDouble()), face.Domain(1).ParameterAt(random.NextDouble()));
                        if (face.IsPointOnFace(point.PositionUV.X, point.PositionUV.Y) > 0)
                        {
                            point.FaceID = index;
                            point.Position = face.PointAt(point.PositionUV.X, point.PositionUV.Y);
                            initialPoints.Add(point);
                        }
                        else j--;
                    }
                }), i as object));
            }
            DoTasksParallel(tasks);




            //Build grid partition
            double radius = distance / 2;
            double cellSize = radius / (Math.Sqrt(3));
            Box boundingBoxPoints = new Box(new BoundingBox(initialPoints.Select(Point => Point.Position)));
            int nX = Convert.ToInt32(boundingBoxPoints.X.Length / cellSize);
            int nY = Convert.ToInt32(boundingBoxPoints.Y.Length / cellSize);
            int nZ = Convert.ToInt32(boundingBoxPoints.Z.Length / cellSize);

            Box boundingBoxCells = new Box(boundingBoxPoints.Plane, new Interval(0, nX * cellSize), new Interval(0, nY * cellSize), new Interval(0, nZ * cellSize));
            boundingBoxCells.Transform(Transform.Translation(boundingBoxPoints.Center - boundingBoxCells.Center));

            int nPhaseGroups = P * P * P;
            int nCellsInPhaseGroups = (nX * nY * nZ) / nPhaseGroups;
            List<ConcurrentBag<Cell>> phaseGroups = new List<ConcurrentBag<Cell>>(nPhaseGroups);
            Cell[,,] cells = new Cell[nX, nY, nZ];

            int nPhaseGroup = 0;
            for (int x = 0; x < P; x++)
            {
                for (int y = 0; y < P; y++)
                {
                    for (int z = 0; z < P; z++)
                    {
                        phaseGroups.Add(new ConcurrentBag<Cell>());
                        for (int px = x; px < nX; px += P)
                        {
                            for (int py = y; py < nY; py += P)
                            {
                                for (int pz = z; pz < nZ; pz += P)
                                {
                                    Cell cell = new Cell();
                                    cell.X = px;
                                    cell.Y = py;
                                    cell.Z = pz;
                                    cell.PhaseGroup = nPhaseGroup;
                                    cells[px, py, pz] = cell;

                                    tasks.Add(new Task(new Action(() =>
                                    {
                                        Point3d corner0 = boundingBoxCells.PointAt(cell.X / (double)nX, cell.Y / (double)nY, cell.Z / (double)nZ);
                                        Point3d corner1 = corner0 + new Point3d(cellSize, cellSize, cellSize);
                                        cell.Box = new Box(boundingBoxCells.Plane, new Point3d[2] { corner0, corner1 });
                                        foreach (var point in initialPoints)
                                        {
                                            if (point.Cell.PhaseGroup > -1) continue;                                            
                                            if(cell.Box.Contains(point.Position))
                                            {
                                                point.Cell = cell;
                                                cell.InitialPoints.Add(point);
                                            }
                                        }
                                        if(cell.InitialPoints.Count > 0)
                                        {
                                            phaseGroups[cell.PhaseGroup].Add(cell);
                                        }
                                    })));
                                }
                            }
                        }
                        nPhaseGroup++;
                    }
                }
            }
            DoTasksParallel(tasks);




            //Parallel sampling
            for (int i = 0; i < K; i++)
            {
                phaseGroups = phaseGroups.OrderBy(a => random.Next()).ToList();
                foreach (var phaseGroup in phaseGroups)
                {
                    foreach (var cell in phaseGroup)
                    {
                        if (i >= cell.InitialPoints.Count) continue;

                        tasks.Add(new Task(new Action(() =>
                        {
                            List<Point> nearPoints = new List<Point>();
                            for (int x = -2; x <= 2; x++)
                            {
                                int iX = cell.X + x;
                                if (iX < 0 || iX >= nX) continue;
                                for (int y = -2; y <= 2; y++)
                                {
                                    int iY = cell.Y + y;
                                    if (iY < 0 || iY >= nY) continue;
                                    for (int z = -2; z <= 2; z++)
                                    {
                                        int iZ = cell.Z + z;
                                        if (iZ < 0 || iZ >= nZ) continue;
                                        nearPoints.AddRange(cells[iX, iY, iZ].Samples);
                                    }
                                }
                            }

                            Point point = cell.InitialPoints[i];
                            point.normal = brep.Faces[point.FaceID].NormalAt(point.PositionUV.X, point.PositionUV.Y);

                            nearPoints.Sort((x, y) =>
                            {
                                return x.Position.DistanceToSquared(point.Position).CompareTo(y.Position.DistanceToSquared(point.Position));
                            });

                            bool valid = true;
                            foreach (var other in nearPoints)
                            {
                                double distanceEuclidien = other.Position.DistanceTo(point.Position);
                                Vector3d vector = (other.Position - point.Position) / distanceEuclidien;
                                double c1 = point.normal * vector;
                                double c2 = other.normal * vector;

                                double distanceGeodesic = ((Math.Asin(c1) - Math.Asin(c2)) / c1 - c2) * distanceEuclidien;

                                if (distanceGeodesic < distance)
                                {
                                    valid = false;
                                    break;
                                }
                            }
                            if (valid) cell.Samples.Add(point);
                        })));
                    }
                    DoTasksParallel(tasks);
                }
            }



            //Collect samples;
            List<Point> samples = new List<Point>();
            foreach (var phaseGroup in phaseGroups)
            {
                foreach (var cell in phaseGroup)
                {
                    samples.AddRange(cell.Samples);
                }
            }

            return samples;
        }



        void DoTasksParallel(List<Task> tasks)
        {
            foreach (var task in tasks)
            {
                task.Start();
            }
            foreach (var task in tasks)
            {
                task.Wait();
            }
            tasks.Clear();
        }
        void DoTasksSequential(List<Task> tasks)
        {
            foreach (var task in tasks)
            {
                task.Start();
                task.Wait();
            }
        }

        /// <summary>
        /// Provides an Icon for the component.
        /// </summary>
        protected override System.Drawing.Bitmap Icon
        {
            get
            {
                //You can add image files to your project resources and access them like this:
                // return Resources.IconForThisComponent;
                return null;
            }
        }

        /// <summary>
        /// Gets the unique ID for this component. Do not change this ID after release.
        /// </summary>
        public override Guid ComponentGuid
        {
            get { return new Guid("8678b2bc-20ed-4539-908e-168440588b02"); }
        }
    }
}