using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Runtime.CompilerServices;
using System.Text;
using System.Threading.Tasks;
using Grasshopper.Kernel;
using Grasshopper.Kernel.Parameters;
using Rhino;
using Rhino.Display;
using Rhino.DocObjects;
using Rhino.Geometry;
using Rhino.Render.UI;

// In order to load the result of this wizard, you will also need to
// add the output bin/ folder of this project to the list of loaded
// folder in Grasshopper.
// You can use the _GrasshopperDeveloperSettings Rhino command for that.

namespace Poisson_Disk_Sampling
{
    public class PoissonDiskSampling3D : GH_Component
    {
        /// <summary>
        /// Each implementation of GH_Component must provide a public 
        /// constructor without any arguments.
        /// Category represents the Tab in which the component will appear, 
        /// Subcategory the panel. If you use non-existing tab or panel names, 
        /// new tabs/panels will automatically be created.
        /// </summary>
        public PoissonDiskSampling3D()
          : base("Poisson Sphere Sampling 3D", "Poisson",
              "Sample points inside a curve with the poisson sphere sampling method.",
              "Vector", "Grid")
        {
        }

        public override GH_Exposure Exposure => GH_Exposure.tertiary | GH_Exposure.obscure;

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            Param.Brep = pManager.AddBrepParameter("Boundary Brep", "B", "The outer boundary to sample points in.", GH_ParamAccess.item);
            Param.Distance = pManager.AddNumberParameter("Distance", "D", "The distance between points or the diameter of a sphere.", GH_ParamAccess.item);
            Param.Seed = pManager.AddIntegerParameter("Seed", "S", "The seed that for random generation of the points.", GH_ParamAccess.item);
            Param.Random = pManager.AddIntegerParameter("Random", "R", "Should the generation be truly random or deterministic?", GH_ParamAccess.item, 0);

            Param_Integer paramRandom = pManager[3] as Param_Integer;
            paramRandom.AddNamedValue("Deterministic", 0);
            paramRandom.AddNamedValue("Random", 1);

            pManager[2].Optional = true;
            pManager[3].Optional = true;
        }

        struct Param
        {
            public static int Brep, Distance, Seed, Random, Samples;
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            Param.Samples = pManager.AddPointParameter("Sampled Points", "P", "Points sampled by the poisson sphere sampling method.", GH_ParamAccess.list);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object can be used to retrieve data from input parameters and 
        /// to store data in output parameters.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            Brep brep = null;
            double distance = 0;
            int seed = 0;
            int random = 0;

            if (!DA.GetData(Param.Brep, ref brep)) return;
            if (!DA.GetData(Param.Distance, ref distance)) return;
            DA.GetData(Param.Seed, ref seed);
            DA.GetData(Param.Random, ref random);
            
            if(random != 0 && random != 1)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, "Please choose a valid selection of randomness.");
            }
            bool cancel = false;
            if(!brep.IsValid)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Brep is not valid!");
                cancel = true;
            }
            if(!brep.IsSolid)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Brep is not a solid!");
                cancel = true;
            }
            if (cancel) return;

            List<Point3d> samples;
            if (random == 1) { samples = SampleRandom(brep, distance, seed); }
            else { samples = Sample(brep, distance, seed); }

            DA.SetDataList(Param.Samples, samples);
        }

        ConcurrentBag<Point3d> samples;
        Brep boundary;
        double distance;
        double distanceSquared;
        double cellMinimum;
        int seed;
        bool random;
        struct SampleCellArguments
        {
            public SampleCellArguments(Box cell, int depth)
            {
                this.cell = cell;
                this.depth = depth;
            }
            public Box cell;
            public int depth;
        }

        public List<Point3d> SampleRandom(Brep brep, double distance, int seed)
        {
            return Sample(brep, distance, seed, true);
        }
        public List<Point3d> Sample(Brep brep, double distance, int seed)
        {
            return Sample(brep, distance, seed, false);
        }

        List<Point3d> Sample(Brep brep, double distance, int seed, bool random)
        {
            if (distance <= 0) return null;

            Plane plane;
            BoundingBox boundingBox = brep.GetBoundingBox(true);

            Box box = new Box(boundingBox);

            //double cellLength = box.X.Length > box.Y.Length ? box.X.Length : box.Y.Length;

            //Box cell = new Box(box.Plane, new Interval(box.X.Min, box.X.Min + cellLength), new Interval(box.Y.Min, box.Y.Min + cellLength), new Interval());
            Box cell = new Box(boundingBox);

            samples = new ConcurrentBag<Point3d>();
            boundary = brep;
            this.distance = distance;
            distanceSquared = distance * distance;
            cellMinimum = (distance/2)/Math.Sqrt(2);
            this.seed = seed;
            this.random = random;

            SampleCell(new SampleCellArguments(cell, 0));

            return new List<Point3d>(samples);
        }

        void SampleCell(object argumentsObject)
        {
            SampleCellArguments arguments = (SampleCellArguments)argumentsObject;
            Box cell = arguments.cell;
            int depth = arguments.depth;

            Random random = new Random(seed + depth);

            Point3d sample = cell.PointAt(random.NextDouble(), random.NextDouble(), random.NextDouble());

            bool valid = true;

            if (sample.DistanceTo(boundary.ClosestPoint(sample)) < distance / 2) valid = false;
            if(valid)
            {
                if (!boundary.IsPointInside(sample, RhinoDoc.ActiveDoc.ModelAbsoluteTolerance, true)) valid = false;
            }
            if (valid)
            {
                foreach (var point in samples)
                {
                    if (sample.DistanceToSquared(point) < distanceSquared)
                    {
                        valid = false;
                        break;
                    }
                }
            }
            if (valid) samples.Add(sample);

            if (cell.X.Length < cellMinimum) return;

            List<Box> newCells = new List<Box>();
            for (double u = 0; u < 1; u += 0.5)
            {
                for (double v = 0; v < 1; v += 0.5)
                {
                    for (double w = 0; w < 1; w += 0.5)
                    {
                        Box newCell = new Box(cell.Plane, new Point3d[] { cell.PointAt(u, v, w), cell.PointAt(u + 0.5, v + 0.5, w + 0.5) });
                        newCells.Add(newCell);
                    }
                }
            }

            List<Task> tasks = new List<Task>(4);
            for (int i = newCells.Count-1; i >= 0; i--)
            {
                Box newCell = newCells[random.Next(i)];
                newCells.Remove(newCell);
                SampleCellArguments newArguments = new SampleCellArguments(newCell, depth + 1);

                if (this.random)
                {
                    Task task = new Task(new Action<object>(SampleCell), newArguments);
                    tasks.Add(task);
                    task.Start();
                }
                else { SampleCell(newArguments); }
            }

            foreach (Task task in tasks)
            {
                task.Wait();
            }
            return;
        }



        /// <summary>
        /// Provides an Icon for every component that will be visible in the User Interface.
        /// Icons need to be 24x24 pixels.
        /// </summary>
        protected override System.Drawing.Bitmap Icon
        {
            get
            {
                // You can add image files to your project resources and access them like this:
                //return Resources.IconForThisComponent;
                return Properties.Resources.PoissonDiskSampling3D;
                //return null;
            }
        }

        /// <summary>
        /// Each component must have a unique Guid to identify it. 
        /// It is vital this Guid doesn't change otherwise old ghx files 
        /// that use the old ID will partially fail during loading.
        /// </summary>
        public override Guid ComponentGuid
        {
            get { return new Guid("631ac774-accd-476a-ba0a-1bac3897d253"); }
        }
    }
}
