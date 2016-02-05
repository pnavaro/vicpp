#include "xdmf.hpp"

Xdmf::Xdmf(const Grid & G)
{
    FILE *xmf = 0;

    int NX = G.nx;
    int NY = G.ny;
    int NZ = G.nz;

    // Open the file and write the XML description of the mesh..
    xmf = fopen("field.xmf", "w");
    fprintf(xmf, "<?xml version=\"1.0\" ?>\n");
    fprintf(xmf, "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n");
    fprintf(xmf, "<Xdmf Version=\"2.0\">\n");
    fprintf(xmf, "<Domain>\n");
    fprintf(xmf, "<Grid Name=\"mesh\" GridType=\"Uniform\">\n");
    fprintf(xmf, "<Topology TopologyType=\"3DSMesh\" NumberOfElements=\"%d %d %d\"/>\n", NX, NY, NZ);
    fprintf(xmf, "<Geometry GeometryType=\"X_Y_Z\">\n");
    fprintf(xmf, "<DataItem Dimensions=\"%d %d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n",NX,NY,NZ);
    fprintf(xmf, "mesh.h5:/X\n");
    fprintf(xmf, "</DataItem>\n");
    fprintf(xmf, "<DataItem Dimensions=\"%d %d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n",NX,NY,NZ);
    fprintf(xmf, "mesh.h5:/Y\n");
    fprintf(xmf, "</DataItem>\n");
    fprintf(xmf, "<DataItem Dimensions=\"%d %d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n",NX,NY,NZ);
    fprintf(xmf, "mesh.h5:/Z\n");
    fprintf(xmf, "</DataItem>\n");
    fprintf(xmf, "</Geometry>\n");
    fprintf(xmf, "<Attribute Name=\"U\" AttributeType=\"Scalar\" Center=\"Node\">\n");
    fprintf(xmf, "<DataItem Dimensions=\"%d %d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n",NX,NY,NZ);
    fprintf(xmf, "field.h5:/U\n");
    fprintf(xmf, "</DataItem>\n");
    fprintf(xmf, "</Attribute>\n");
    fprintf(xmf, "<Attribute Name=\"V\" AttributeType=\"Scalar\" Center=\"Node\">\n");
    fprintf(xmf, "<DataItem Dimensions=\"%d %d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n",NX,NY,NZ);
    fprintf(xmf, "field.h5:/V\n");
    fprintf(xmf, "</DataItem>\n");
    fprintf(xmf, "</Attribute>\n");
    fprintf(xmf, "</Grid>\n");
    fprintf(xmf, "</Domain>\n");
    fprintf(xmf, "</Xdmf>\n");
    fclose(xmf);


}

void Xdmf::write_hdf5_mesh( const Grid &G)
{
    hid_t       file, dataset;    /* file and dataset handles */
    hid_t       dataspace;        /* handles                  */
    hsize_t     dimsf[3];         /* dataset dimensions       */
    herr_t      status;

    file = H5Fcreate("mesh.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    dimsf[0] = G.nx;
    dimsf[1] = G.ny;
    dimsf[2] = G.ny;

    dataspace = H5Screate_simple(3, dimsf, NULL);

    dataset = H5Dcreate2(file, "/X", H5T_NATIVE_DOUBLE, dataspace,
                         H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &G.x[0]);
    H5Dclose(dataset);

    dataset = H5Dcreate2(file, "/Y", H5T_NATIVE_DOUBLE, dataspace,
                         H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &G.y[0]);
    H5Dclose(dataset);

    dataset = H5Dcreate2(file, "/Z", H5T_NATIVE_DOUBLE, dataspace,
                         H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &G.z[0]);
    H5Dclose(dataset);

    H5Sclose(dataspace);

    H5Fclose(file);

}




void Xdmf::write_hdf5_data(const Grid &G, double * u, double * v)
{

    hid_t       file, dataset;    /* file and dataset handles */
    hid_t       dataspace;        /* handles                  */
    hsize_t     dimsf[3];         /* dataset dimensions       */
    herr_t      status;

    file = H5Fcreate(H5FILE_NAME, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    dimsf[0] = G.nx;
    dimsf[1] = G.ny;
    dimsf[2] = G.nz;

    dataspace = H5Screate_simple(3, dimsf, NULL);
    dataset = H5Dcreate2(file, "/U", H5T_NATIVE_DOUBLE, dataspace,
                         H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &u[0]);
    H5Dclose(dataset);

    dataset = H5Dcreate2(file, "/V", H5T_NATIVE_DOUBLE, dataspace,
                         H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &v[0]);
    H5Dclose(dataset);

    H5Fclose(file);

}
