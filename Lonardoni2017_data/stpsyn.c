/* Created by Language version: 7.7.0 */
/* NOT VECTORIZED */
#define NRN_VECTORIZED 0
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "scoplib_ansi.h"
#undef PI
#define nil 0
#include "md1redef.h"
#include "section.h"
#include "nrniv_mf.h"
#include "md2redef.h"
 
#if METHOD3
extern int _method3;
#endif

#if !NRNGPU
#undef exp
#define exp hoc_Exp
extern double hoc_Exp(double);
#endif
 
#define nrn_init _nrn_init__stp
#define _nrn_initial _nrn_initial__stp
#define nrn_cur _nrn_cur__stp
#define _nrn_current _nrn_current__stp
#define nrn_jacob _nrn_jacob__stp
#define nrn_state _nrn_state__stp
#define _net_receive _net_receive__stp 
 
#define _threadargscomma_ /**/
#define _threadargsprotocomma_ /**/
#define _threadargs_ /**/
#define _threadargsproto_ /**/
 	/*SUPPRESS 761*/
	/*SUPPRESS 762*/
	/*SUPPRESS 763*/
	/*SUPPRESS 765*/
	 extern double *getarg();
 static double *_p; static Datum *_ppvar;
 
#define t nrn_threads->_t
#define dt nrn_threads->_dt
#define tau_1 _p[0]
#define tau_rec _p[1]
#define tau_facil _p[2]
#define U _p[3]
#define u0 _p[4]
#define flag _p[5]
#define z0 _p[6]
#define y0 _p[7]
#define x _p[8]
#define _tsav _p[9]
#define _nd_area  *_ppvar[0]._pval
#define new	*_ppvar[2]._pval
#define _p_new	_ppvar[2]._pval
 
#if MAC
#if !defined(v)
#define v _mlhv
#endif
#if !defined(h)
#define h _mlhh
#endif
#endif
 
#if defined(__cplusplus)
extern "C" {
#endif
 static int hoc_nrnpointerindex =  2;
 /* external NEURON variables */
 /* declaration of user functions */
 static int _mechtype;
extern void _nrn_cacheloop_reg(int, int);
extern void hoc_register_prop_size(int, int, int);
extern void hoc_register_limits(int, HocParmLimits*);
extern void hoc_register_units(int, HocParmUnits*);
extern void nrn_promote(Prop*, int, int);
extern Memb_func* memb_func;
 
#define NMODL_TEXT 1
#if NMODL_TEXT
static const char* nmodl_file_text;
static const char* nmodl_filename;
extern void hoc_reg_nmodl_text(int, const char*);
extern void hoc_reg_nmodl_filename(int, const char*);
#endif

 extern Prop* nrn_point_prop_;
 static int _pointtype;
 static void* _hoc_create_pnt(_ho) Object* _ho; { void* create_point_process();
 return create_point_process(_pointtype, _ho);
}
 static void _hoc_destroy_pnt();
 static double _hoc_loc_pnt(_vptr) void* _vptr; {double loc_point_process();
 return loc_point_process(_pointtype, _vptr);
}
 static double _hoc_has_loc(_vptr) void* _vptr; {double has_loc_point();
 return has_loc_point(_vptr);
}
 static double _hoc_get_loc_pnt(_vptr)void* _vptr; {
 double get_loc_point_process(); return (get_loc_point_process(_vptr));
}
 extern void _nrn_setdata_reg(int, void(*)(Prop*));
 static void _setdata(Prop* _prop) {
 _p = _prop->param; _ppvar = _prop->dparam;
 }
 static void _hoc_setdata(void* _vptr) { Prop* _prop;
 _prop = ((Point_process*)_vptr)->_prop;
   _setdata(_prop);
 }
 /* connect user functions to hoc names */
 static VoidFunc hoc_intfunc[] = {
 0,0
};
 static Member_func _member_func[] = {
 "loc", _hoc_loc_pnt,
 "has_loc", _hoc_has_loc,
 "get_loc", _hoc_get_loc_pnt,
 0, 0
};
 /* declare global and static user variables */
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 "U", 0, 1,
 "tau_facil", 0, 1e+009,
 "tau_rec", 1e-009, 1e+009,
 "tau_1", 1e-009, 1e+009,
 "u0", 0, 1,
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "tau_1", "ms",
 "tau_rec", "ms",
 "tau_facil", "ms",
 "U", "1",
 "u0", "1",
 0,0
};
 static double v = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 0,0
};
 static DoubVec hoc_vdoub[] = {
 0,0,0
};
 static double _sav_indep;
 static void nrn_alloc(Prop*);
static void  nrn_init(_NrnThread*, _Memb_list*, int);
static void nrn_state(_NrnThread*, _Memb_list*, int);
 static void _hoc_destroy_pnt(_vptr) void* _vptr; {
   destroy_point_process(_vptr);
}
 /* connect range variables in _p that hoc is supposed to know about */
 static const char *_mechanism[] = {
 "7.7.0",
"stp",
 "tau_1",
 "tau_rec",
 "tau_facil",
 "U",
 "u0",
 "flag",
 "z0",
 "y0",
 0,
 0,
 0,
 "new",
 0};
 
extern Prop* need_memb(Symbol*);

static void nrn_alloc(Prop* _prop) {
	Prop *prop_ion;
	double *_p; Datum *_ppvar;
  if (nrn_point_prop_) {
	_prop->_alloc_seq = nrn_point_prop_->_alloc_seq;
	_p = nrn_point_prop_->param;
	_ppvar = nrn_point_prop_->dparam;
 }else{
 	_p = nrn_prop_data_alloc(_mechtype, 10, _prop);
 	/*initialize range parameters*/
 	tau_1 = 3;
 	tau_rec = 500;
 	tau_facil = 100;
 	U = 0.5;
 	u0 = 0;
 	flag = 0;
 	z0 = 0.4;
 	y0 = 0.4;
  }
 	_prop->param = _p;
 	_prop->param_size = 10;
  if (!nrn_point_prop_) {
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 3, _prop);
  }
 	_prop->dparam = _ppvar;
 	/*connect ionic variables to this model*/
 
}
 static void _initlists();
 static void _net_receive(Point_process*, double*, double);
 static void _net_init(Point_process*, double*, double);
 extern Symbol* hoc_lookup(const char*);
extern void _nrn_thread_reg(int, int, void(*)(Datum*));
extern void _nrn_thread_table_reg(int, void(*)(double*, Datum*, Datum*, _NrnThread*, int));
extern void hoc_register_tolerance(int, HocStateTolerance*, Symbol***);
extern void _cvode_abstol( Symbol**, double*, int);

 void _stpsyn_reg() {
	int _vectorized = 0;
  _initlists();
 	_pointtype = point_register_mech(_mechanism,
	 nrn_alloc,(void*)0, (void*)0, (void*)0, nrn_init,
	 hoc_nrnpointerindex, 0,
	 _hoc_create_pnt, _hoc_destroy_pnt, _member_func);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_setdata_reg(_mechtype, _setdata);
 #if NMODL_TEXT
  hoc_reg_nmodl_text(_mechtype, nmodl_file_text);
  hoc_reg_nmodl_filename(_mechtype, nmodl_filename);
#endif
  hoc_register_prop_size(_mechtype, 10, 3);
  hoc_register_dparam_semantics(_mechtype, 0, "area");
  hoc_register_dparam_semantics(_mechtype, 1, "pntproc");
  hoc_register_dparam_semantics(_mechtype, 2, "pointer");
 add_nrn_has_net_event(_mechtype);
 pnt_receive[_mechtype] = _net_receive;
 pnt_receive_init[_mechtype] = _net_init;
 pnt_receive_size[_mechtype] = 5;
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 stp C:/Users/teppola/Documents/OCNC_Project/Heidis_modeling_project_OCNC/data-analysis/Lonardoni2017_data/stpsyn.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
static int _reset;
static char *modelname = "";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
 
/*VERBATIM*/
#include <stdlib.h> //open the standard library of C
 
static void _net_receive (_pnt, _args, _lflag) Point_process* _pnt; double* _args; double _lflag; 
{    _p = _pnt->_prop->param; _ppvar = _pnt->_prop->dparam;
  if (_tsav > t){ extern char* hoc_object_name(); hoc_execerror(hoc_object_name(_pnt->ob), ":Event arrived out of order. Must call ParallelContext.set_maxstep AFTER assigning minimum NetCon.delay");}
 _tsav = t; {
   _args[2] = _args[2] * exp ( - ( t - _args[4] ) / tau_rec ) ;
   _args[2] = _args[2] + ( _args[1] * ( exp ( - ( t - _args[4] ) / tau_1 ) - exp ( - ( t - _args[4] ) / tau_rec ) ) / ( ( tau_1 / tau_rec ) - 1.0 ) ) ;
   _args[1] = _args[1] * exp ( - ( t - _args[4] ) / tau_1 ) ;
   x = 1.0 - _args[1] - _args[2] ;
   if ( tau_facil > 0.0 ) {
     _args[3] = u0 + ( _args[3] - u0 ) * exp ( - ( t - _args[4] ) / tau_facil ) ;
     }
   else {
     _args[3] = U ;
     }
   if ( tau_facil > 0.0 ) {
     _args[3] = _args[3] + U * ( 1.0 - _args[3] ) ;
     }
   new = _args[0] * ( x * _args[3] ) ;
   _args[1] = _args[1] + x * _args[3] ;
   net_event ( _pnt, t ) ;
   if ( _lflag > 0.0 ) {
     printf ( "%g -> 0: \t t=%.0g\t new=%.4g\t y=%.4g\t xu=%.4g\t x=%.4g\n" , _lflag , t , new , _args[1] , x * _args[3] , x ) ;
     }
   _args[4] = t ;
   } }
 
static void _net_init(Point_process* _pnt, double* _args, double _lflag) {
    _args[1] = y0 ;
   _args[2] = z0 ;
   _args[3] = u0 ;
   _args[4] = 1000.0 ;
   }

static void initmodel() {
  int _i; double _save;_ninits++;
{

}
}

static void nrn_init(_NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; double _v; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
 _tsav = -1e20;
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 v = _v;
 initmodel();
}}

static double _nrn_current(double _v){double _current=0.;v=_v;{
} return _current;
}

static void nrn_state(_NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; double _v = 0.0; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
 _nd = _ml->_nodelist[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 v=_v;
{
}}

}

static void terminal(){}

static void _initlists() {
 int _i; static int _first = 1;
  if (!_first) return;
_first = 0;
}

#if NMODL_TEXT
static const char* nmodl_filename = "stpsyn.mod";
static const char* nmodl_file_text = 
  "COMMENT\n"
  "Revised 12/15/2000 in light of a personal communication\n"
  "from Misha Tsodyks that u is incremented _before_ x is\n"
  "converted to y--a point that was not clear in the paper.\n"
  "If u is incremented _after_ x is converted to y, then\n"
  "the first synaptic activation after a long interval of\n"
  "silence will produce smaller and smaller postsynaptic\n"
  "effect as the length of the silent interval increases,\n"
  "eventually becoming vanishingly small.\n"
  "\n"
  "Implementation of a model of short-term facilitation and depression\n"
  "based on the kinetics described in\n"
  "  Tsodyks et al.\n"
  "  Synchrony generation in recurrent networks\n"
  "  with frequency-dependent synapses\n"
  "  Journal of Neuroscience 20:RC50:1-5, 2000.\n"
  "Their mechanism represented synapses as current sources.\n"
  "The mechanism implemented here uses a conductance change instead.\n"
  "\n"
  "The basic scheme is\n"
  "\n"
  "x -------> y    Instantaneous, spike triggered.\n"
  "                Increment is u*x (see discussion of u below).\n"
  "                x == fraction of \"synaptic resources\" that have\n"
  "                     \"recovered\" (fraction of xmtr pool that is\n"
  "                     ready for release, or fraction of postsynaptic\n"
  "                     channels that are ready to be opened, or some\n"
  "                     joint function of these two factors)\n"
  "                y == fraction of \"synaptic resources\" that are in the\n"
  "                     \"active state.\"  This is proportional to the\n"
  "                     number of channels that are open, or the\n"
  "                     fraction of max synaptic current that is\n"
  "                     being delivered.\n"
  "  tau_1\n"
  "y -------> z    z == fraction of \"synaptic resources\" that are\n"
  "                     in the \"inactive state\"\n"
  "\n"
  "  tau_rec\n"
  "z -------> x\n"
  "\n"
  "where x + y + z = 1\n"
  "\n"
  "The active state y is multiplied by a synaptic weight to compute\n"
  "the actual synaptic conductance (or current, in the original form\n"
  "of the model).\n"
  "\n"
  "In addition, there is a \"facilition\" term u that\n"
  "governs the fraction of x that is converted to y\n"
  "on each synaptic activation.\n"
  "\n"
  "  -------> u    Instantaneous, spike triggered,\n"
  "                happens _BEFORE_ x is converted to y.\n"
  "                Increment is U*(1-u) where U and u both\n"
  "                lie in the range 0 - 1.\n"
  "  tau_facil\n"
  "u ------->      decay of facilitation\n"
  "\n"
  "This implementation for NEURON offers the user a parameter\n"
  "u0 that has a default value of 0 but can be used to specify\n"
  "a nonzero initial value for u.\n"
  "\n"
  "When tau_facil = 0, u is supposed to equal U.\n"
  "\n"
  "Note that the synaptic conductance in this mechanism\n"
  "has the same kinetics as y, i.e. decays with time\n"
  "constant tau_1.\n"
  "\n"
  "This mechanism can receive multiple streams of\n"
  "synaptic input via NetCon objects.\n"
  "Each stream keeps track of its own\n"
  "weight and activation history.\n"
  "\n"
  "The printf() statements are for testing purposes only.\n"
  "ENDCOMMENT\n"
  "\n"
  "VERBATIM\n"
  "#include <stdlib.h> //open the standard library of C\n"
  "ENDVERBATIM\n"
  "\n"
  "NEURON {\n"
  "	POINT_PROCESS stp\n"
  "    POINTER new\n"
  "	RANGE tau_1, tau_rec, tau_facil, U, u0,flag,z0,y0\n"
  "\n"
  "}\n"
  "\n"
  "UNITS {\n"
  "	(nA) = (nanoamp)\n"
  "	(mV) = (millivolt)\n"
  "	(umho) = (micromho)\n"
  "}\n"
  "\n"
  "PARAMETER {\n"
  "	: e = -90 mV for inhibitory synapses, 0 mV for excitatory\n"
  "	: tau_1 was the same for inhibitory and excitatory synapses\n"
  "	: in the models used by T et al.\n"
  "  : tau_1 = 3 ms\n"
  "  : tau_1 = 1 ms after fitting\n"
  "	: tau_rec = 100 ms for inhibitory synapses,\n"
  "	:           800 ms for excitatory\n"
  "	: tau_facil = 100 ms for inhibitory synapses,\n"
  "	:             0 ms for excitatory\n"
  "  : tau_facil = 200 after fitting\n"
  "	: U = 0.04 for inhibitory synapses, 0.5 for excitatory\n"
  "	: the (1) is needed for the < 0, 1 > to be effective\n"
  "	:   in limiting the values of U and u0\n"
  "  : 0.15 in the ctrl case, 0.36 in the Abeta case\n"
  "	: initial value for the \"facilitation variable\"\n"
  "  : u0 = 0 Abeta Tsodyks 1.1\n"
  "\n"
  "\n"
  "	tau_1 = 3 (ms) < 1e-9, 1e9 >\n"
  "	tau_rec = 500 (ms) < 1e-9, 1e9 >\n"
  "	tau_facil = 100 (ms) < 0, 1e9 >\n"
  "	U = 0.5 (1) < 0, 1 >\n"
  "	u0 = 0 (1) < 0, 1 >\n"
  "	flag=0\n"
  "	z0=.4\n"
  "	y0=.4\n"
  "}\n"
  "\n"
  "ASSIGNED {\n"
  "	x\n"
  "	new\n"
  "}\n"
  "\n"
  "NET_RECEIVE(weight (umho), y, z, u, tsyn (ms)) {\n"
  "  INITIAL {\n"
  "  : these are in NET_RECEIVE to be per-stream\n"
  "  	y = y0\n"
  "  	z = z0\n"
  "  :	u = 0\n"
  "  	u = u0\n"
  "  	tsyn = 1000\n"
  "  : this header will appear once per stream\n"
  "   :printf(\"t\\t t-tsyn\\t y\\t z\\t u\\t newu\\t g\\t dg\\t newg\\t newy\\n\")\n"
  "  }\n"
  "\n"
  "	: first calculate z at event-\n"
  "	: based on prior y and z\n"
  "	z = z*exp(-(t - tsyn)/tau_rec)\n"
  "	z = z + ( y*(exp(-(t - tsyn)/tau_1) - exp(-(t - tsyn)/tau_rec)) / ((tau_1/tau_rec)-1) )\n"
  "	: now calc y at event-\n"
  "	y = y*exp(-(t - tsyn)/tau_1)\n"
  "\n"
  "	x = 1-y-z\n"
  "\n"
  "	: calc u at event--\n"
  "	if (tau_facil > 0) {\n"
  "		u = u0 + (u-u0)*exp(-(t - tsyn)/tau_facil)\n"
  "	} else {\n"
  "		u = U\n"
  "	}\n"
  "\n"
  "\n"
  "\n"
  "	if (tau_facil > 0) {\n"
  "		u= u + U*(1-u)\n"
  "	}\n"
  "    new=weight*(x*u)\n"
  "	y=y + x*u\n"
  "    net_event(t)\n"
  "    if (flag>0){\n"
  "    printf(\"%g -> 0: \\t t=%.0g\\t new=%.4g\\t y=%.4g\\t xu=%.4g\\t x=%.4g\\n\",flag, t, new, y, x*u, x)}\n"
  "	tsyn = t\n"
  "\n"
  "\n"
  "}\n"
  ;
#endif
