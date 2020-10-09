#include "Python.h"

extern void shock_csm(char*);
extern void rad_transfer_csm(double, char*, char*);

// definition of shock flux calculator method
static PyObject* lightcurve_shock(PyObject* self, PyObject* args, PyObject* kw)
{
	const char* file_input = NULL;
	static char* argnames[] = {"file_input", NULL};
	if (!PyArg_ParseTupleAndKeywords(args, kw, "|s", argnames, &file_input))
		return NULL;
	shock_csm(file_input);
	return Py_BuildValue("");
}

// definition of transfer method
static PyObject* lightcurve_transfer(PyObject* self, PyObject* args, PyObject* kw)
{
	double r_out;
	const char* file_csm=NULL;
	const char* file_shock=NULL;
	static char* argnames[] = {"r_out", "file_csm", "file_shock", NULL};
	if (!PyArg_ParseTupleAndKeywords(args, kw, "d|ss", argnames, &r_out, &file_csm, &file_shock))
		return NULL;
	rad_transfer_csm(r_out, file_csm, file_shock);
	return Py_BuildValue("");
}

// definition of all methods of this module
static PyMethodDef lightcurve_methods[] = {
	{"shock", lightcurve_shock, METH_VARARGS | METH_KEYWORDS},
	{"transfer", lightcurve_transfer, METH_VARARGS | METH_KEYWORDS},
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
