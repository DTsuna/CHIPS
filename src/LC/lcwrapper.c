#include "Python.h"

extern void shock_csm(double, double, double, double, const char*, const char*);
extern void rad_transfer_csm(double, double, double, double, double, const char*, const char*, const char*, const char*, const char*);

// definition of shock flux calculator method
static PyObject* lightcurve_shock(PyObject* self, PyObject* args, PyObject* kw)
{
	double E_ej, M_ej, n, delta;
	const char* file_csm = NULL;
	const char* file_output = NULL;
	static char* argnames[] = {"E_ej", "M_ej", "n", "delta", "file_csm", "file_output", NULL};
	if (!PyArg_ParseTupleAndKeywords(args, kw, "dddd|ss", argnames, &E_ej, &M_ej, &n, &delta, &file_csm, &file_output))
		return NULL;
	shock_csm(E_ej, M_ej, n, delta, file_csm, file_output);
	return Py_BuildValue("");
}

// definition of transfer method
static PyObject* lightcurve_transfer(PyObject* self, PyObject* args, PyObject* kw)
{
	double r_out;
	double Eej, Mej, n, delta;
	const char* file_csm=NULL;
	const char* file_shock=NULL;
	const char* file_lc=NULL;
	const char* file_lc_band=NULL;
	const char* dir_name_shockprofiles = NULL;
	static char* argnames[] = {"Eej", "Mej", "n", "delta", "r_out", "file_csm", "file_shock", "file_lc", "file_lc_band", "dir_name_shockprofiles", NULL};
	if (!PyArg_ParseTupleAndKeywords(args, kw, "ddddd|sssss", argnames, &Eej, &Mej, &n, &delta, &r_out, &file_csm, &file_shock, &file_lc, &file_lc_band, &dir_name_shockprofiles))
		return NULL;
	rad_transfer_csm(Eej, Mej, n, delta, r_out, file_csm, file_shock, file_lc, file_lc_band, dir_name_shockprofiles);
	return Py_BuildValue("");
}

// definition of all methods of this module
static PyMethodDef lightcurve_methods[] = {
	{"shock", (PyCFunction) lightcurve_shock, METH_VARARGS | METH_KEYWORDS},
	{"transfer", (PyCFunction) lightcurve_transfer, METH_VARARGS | METH_KEYWORDS},
	{NULL, NULL}
};

// module definition struct
static struct PyModuleDef lightcurve = {
	PyModuleDef_HEAD_INIT,
	"lightcurve",
	"Python3 C API Module(Sample 1)",
	-1,
	lightcurve_methods
};

PyMODINIT_FUNC PyInit_lightcurve(void)
{
	return PyModule_Create(&lightcurve);
}
