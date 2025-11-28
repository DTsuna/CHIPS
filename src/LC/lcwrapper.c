#include "Python.h"

extern void shock_csm(double, double, double, double, double, double, const char*, const char*, int);
extern void rad_transfer_csm(double, double, double, double, double, double, double, 
				const char*, const char*, const char*, const char*, const char*, int);
extern void gen_muTem_table(int);

// definition of shock flux calculator method
static PyObject* lightcurve_shock(PyObject* self, PyObject* args, PyObject* kw)
{
	double E_ej, M_ej, M_ni, n, delta, s;
	int D;
	const char* file_csm = NULL;
	const char* file_output = NULL;
	static char* argnames[] = {"E_ej", "M_ej", "M_ni", "n", "delta", "s", "file_csm", "file_output", "D", NULL};
	if (!PyArg_ParseTupleAndKeywords(args, kw, "dddddd|ssi", argnames, &E_ej, &M_ej, &M_ni, &n, &delta, &s, &file_csm, &file_output, &D))
		return NULL;
	shock_csm(E_ej, M_ej, M_ni, n, delta, s, file_csm, file_output, D);
	return Py_BuildValue("");
}

// generate tables for mu, T. This might be deleted after debug.
static PyObject* lightcurve_opacTable(PyObject* self, PyObject* args, PyObject* kw)
{
	int discriminant;
	static char* argnames[] = {"discriminant", NULL};
	if (!PyArg_ParseTupleAndKeywords(args, kw, "i", argnames, &discriminant))
		return NULL;
	gen_muTem_table(discriminant);
	return Py_BuildValue("");
}

// definition of transfer method
static PyObject* lightcurve_transfer(PyObject* self, PyObject* args, PyObject* kw)
{
	double r_out;
	double Eej, Mej, Mni, n, delta, s;
	int D;
	const char* file_csm=NULL;
	const char* file_shock=NULL;
	const char* file_lc=NULL;
	const char* file_lc_band=NULL;
	const char* dir_Lnu = NULL;
	static char* argnames[] = {"Eej", "Mej", "Mni", "n", "delta", "r_out", "s", "file_csm", "file_shock", "file_lc", "file_lc_band", "dir_Lnu", "D", NULL};
	if (!PyArg_ParseTupleAndKeywords(args, kw, "ddddddd|sssssi", argnames, &Eej, &Mej, &Mni, &n, &delta, &r_out, &s, &file_csm, &file_shock, &file_lc, &file_lc_band, &dir_Lnu, &D))
		return NULL;
	rad_transfer_csm(Eej, Mej, Mni, n, delta, r_out, s, file_csm, file_shock, file_lc, file_lc_band, dir_Lnu, D);
	return Py_BuildValue("");
}

// definition of all methods of this module
static PyMethodDef lightcurve_methods[] = {
	{"shock", (PyCFunction) lightcurve_shock, METH_VARARGS | METH_KEYWORDS},
	{"transfer", (PyCFunction) lightcurve_transfer, METH_VARARGS | METH_KEYWORDS},
	{"opacTable", (PyCFunction) lightcurve_opacTable, METH_VARARGS | METH_KEYWORDS},
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
