using System;
using System.Drawing;
using Grasshopper.Kernel;

namespace Poisson_Disk_Sampling
{
    public class PoissonDiskSamplingInfo : GH_AssemblyInfo
  {
    public override string Name
    {
        get
        {
            return "PoissonDiskSampling";
        }
    }
    public override Bitmap Icon
    {
        get
        {
            //Return a 24x24 pixel bitmap to represent this GHA library.
            return null;
        }
    }
    public override string Description
    {
        get
        {
            //Return a short string describing the purpose of this GHA library.
            return "Sample points with the poisson disk sampling method.";
        }
    }
    public override Guid Id
    {
        get
        {
            return new Guid("34a0e77a-925a-4d8e-aa54-3d37cd7b91b2");
        }
    }

    public override string AuthorName
    {
        get
        {
            //Return a string identifying you or your company.
            return "Peter Krattenmacher";
        }
    }
    public override string AuthorContact
    {
        get
        {
            //Return a string representing your preferred contact details.
            return "mail+grasshoppper@pkratten.de";
        }
    }
}
}
