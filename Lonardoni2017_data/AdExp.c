/* Created by Language version: 7.7.0 */
/* VECTORIZED */
#define NRN_VECTORIZED 1
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
 
#define nrn_init _nrn_init__AdExp
#define _nrn_initial _nrn_initial__AdExp
#define nrn_cur _nrn_cur__AdExp
#define _nrn_current _nrn_current__AdExp
#define nrn_jacob _nrn_jacob__AdExp
#define nrn_state _nrn_state__AdExp
#define _net_receive _net_receive__AdExp 
#define states states__AdExp 
 
#define _threadargscomma_ _p, _ppvar, _thread, _nt,
#define _threadargsprotocomma_ double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt,
#define _threadargs_ _p, _ppvar, _thread, _nt
#define _threadargsproto_ double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt
 	/*SUPPRESS 761*/
	/*SUPPRESS 762*/
	/*SUPPRESS 763*/
	/*SUPPRESS 765*/
	 extern double *getarg();
 /* Thread safe. No static _p or _ppvar. */
 
#define t _nt->_t
#define dt _nt->_dt
#define num _p[0]
#define V_reset _p[1]
#define V_thre _p[2]
#define a _p[3]
#define b _p[4]
#define tau_w _p[5]
#define E_l _p[6]
#define G_l _p[7]
#define Delta_T _p[8]
#define iEXT _p[9]
#define C _p[10]
#define tauEXC _p[11]
#define tauINH _p[12]
#define erev _p[13]
#define irev _p[14]
#define tRISE _p[15]
#define mNMDA _p[16]
#define mAMPA _p[17]
#define v0_block _p[18]
#define k_block _p[19]
#define tRISEnmda _p[20]
#define tDECAYnmda _p[21]
#define Erev _p[22]
#define rpeso _p[23]
#define gNMDA _p[24]
#define gAMPA _p[25]
#define gGABA _p[26]
#define iAMPA _p[27]
#define iGABA _p[28]
#define iNMDA _p[29]
#define iTotal _p[30]
#define vv _p[31]
#define ww _p[32]
#define gEXC _p[33]
#define gINH _p[34]
#define gRISEexc _p[35]
#define gRISEinh _p[36]
#define y1 _p[37]
#define y2 _p[38]
#define iNOISE _p[39]
#define Dvv _p[40]
#define Dww _p[41]
#define DgEXC _p[42]
#define DgINH _p[43]
#define DgRISEexc _p[44]
#define DgRISEinh _p[45]
#define Dy1 _p[46]
#define Dy2 _p[47]
#define v _p[48]
#define _g _p[49]
#define _tsav _p[50]
#define _nd_area  *_ppvar[0]._pval
 
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
 static int hoc_nrnpointerindex =  -1;
 static Datum* _extcall_thread;
 static Prop* _extcall_prop;
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
 _extcall_prop = _prop;
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
#define V_spike V_spike_AdExp
 double V_spike = -10;
#define eps eps_AdExp
 double eps = 0.1;
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "V_spike_AdExp", "mV",
 "eps_AdExp", "ms",
 "V_reset", "mV",
 "V_thre", "mV",
 "a", "nS",
 "b", "pA",
 "tau_w", "ms",
 "E_l", "mV",
 "G_l", "nS",
 "Delta_T", "mV",
 "iEXT", "pA",
 "C", "nS",
 "tauEXC", "ms",
 "tauINH", "ms",
 "erev", "mV",
 "irev", "mV",
 "tRISE", "ms",
 "tRISEnmda", "ms",
 "tDECAYnmda", "ms",
 "Erev", "mV",
 "rpeso", "nS",
 "vv", "mV",
 "ww", "pA",
 "gEXC", "nS",
 "gINH", "nS",
 "gRISEexc", "nS",
 "gRISEinh", "nS",
 "y1", "nS",
 "y2", "nS",
 "gNMDA", "nS",
 "gAMPA", "nS",
 "gGABA", "nS",
 "iAMPA", "pA",
 "iGABA", "pA",
 "iNMDA", "pA",
 "iTotal", "pA",
 0,0
};
 static double delta_t = 0.01;
 static double gRISEinh0 = 0;
 static double gRISEexc0 = 0;
 static double gINH0 = 0;
 static double gEXC0 = 0;
 static double vv0 = 0;
 static double ww0 = 0;
 static double y20 = 0;
 static double y10 = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 "V_spike_AdExp", &V_spike_AdExp,
 "eps_AdExp", &eps_AdExp,
 0,0
};
 static DoubVec hoc_vdoub[] = {
 0,0,0
};
 static double _sav_indep;
 static void nrn_alloc(Prop*);
static void  nrn_init(_NrnThread*, _Memb_list*, int);
static void nrn_state(_NrnThread*, _Memb_list*, int);
 static void nrn_cur(_NrnThread*, _Memb_list*, int);
static void  nrn_jacob(_NrnThread*, _Memb_list*, int);
 
#define _watch_array _ppvar + 3 
 static void _hoc_destroy_pnt(_vptr) void* _vptr; {
   Prop* _prop = ((Point_process*)_vptr)->_prop;
   if (_prop) { _nrn_free_watch(_prop->dparam, 3, 2);}
   destroy_point_process(_vptr);
}
 
static int _ode_count(int);
static void _ode_map(int, double**, double**, double*, Datum*, double*, int);
static void _ode_spec(_NrnThread*, _Memb_list*, int);
static void _ode_matsol(_NrnThread*, _Memb_list*, int);
 
#define _cvode_ieq _ppvar[5]._i
 static void _ode_matsol_instance1(_threadargsproto_);
 /* connect range variables in _p that hoc is supposed to know about */
 static const char *_mechanism[] = {
 "7.7.0",
"AdExp",
 "num",
 "V_reset",
 "V_thre",
 "a",
 "b",
 "tau_w",
 "E_l",
 "G_l",
 "Delta_T",
 "iEXT",
 "C",
 "tauEXC",
 "tauINH",
 "erev",
 "irev",
 "tRISE",
 "mNMDA",
 "mAMPA",
 "v0_block",
 "k_block",
 "tRISEnmda",
 "tDECAYnmda",
 "Erev",
 "rpeso",
 0,
 "gNMDA",
 "gAMPA",
 "gGABA",
 "iAMPA",
 "iGABA",
 "iNMDA",
 "iTotal",
 0,
 "vv",
 "ww",
 "gEXC",
 "gINH",
 "gRISEexc",
 "gRISEinh",
 "y1",
 "y2",
 0,
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
 	_p = nrn_prop_data_alloc(_mechtype, 51, _prop);
 	/*initialize range parameters*/
 	num = -1;
 	V_reset = -60;
 	V_thre = -50.4;
 	a = 4;
 	b = 80.5;
 	tau_w = 144;
 	E_l = -70.6;
 	G_l = 30;
 	Delta_T = 2;
 	iEXT = 0;
 	C = 281;
 	tauEXC = 3;
 	tauINH = 8;
 	erev = 0;
 	irev = -70;
 	tRISE = 1;
 	mNMDA = 0.05;
 	mAMPA = 1;
 	v0_block = -50;
 	k_block = 8;
 	tRISEnmda = 5.63;
 	tDECAYnmda = 140;
 	Erev = 0;
 	rpeso = 12;
  }
 	_prop->param = _p;
 	_prop->param_size = 51;
  if (!nrn_point_prop_) {
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 6, _prop);
  }
 	_prop->dparam = _ppvar;
 	/*connect ionic variables to this model*/
 
}
 static void _initlists();
  /* some states have an absolute tolerance */
 static Symbol** _atollist;
 static HocStateTolerance _hoc_state_tol[] = {
 0,0
};
 
#define _tqitem &(_ppvar[2]._pvoid)
 static void _net_receive(Point_process*, double*, double);
 static void _net_init(Point_process*, double*, double);
 static void _thread_mem_init(Datum*);
 static void _thread_cleanup(Datum*);
 extern Symbol* hoc_lookup(const char*);
extern void _nrn_thread_reg(int, int, void(*)(Datum*));
extern void _nrn_thread_table_reg(int, void(*)(double*, Datum*, Datum*, _NrnThread*, int));
extern void hoc_register_tolerance(int, HocStateTolerance*, Symbol***);
extern void _cvode_abstol( Symbol**, double*, int);

 void _AdExp_reg() {
	int _vectorized = 1;
  _initlists();
 	_pointtype = point_register_mech(_mechanism,
	 nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init,
	 hoc_nrnpointerindex, 5,
	 _hoc_create_pnt, _hoc_destroy_pnt, _member_func);
  _extcall_thread = (Datum*)ecalloc(4, sizeof(Datum));
  _thread_mem_init(_extcall_thread);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_setdata_reg(_mechtype, _setdata);
     _nrn_thread_reg(_mechtype, 1, _thread_mem_init);
     _nrn_thread_reg(_mechtype, 0, _thread_cleanup);
 #if NMODL_TEXT
  hoc_reg_nmodl_text(_mechtype, nmodl_file_text);
  hoc_reg_nmodl_filename(_mechtype, nmodl_filename);
#endif
  hoc_register_prop_size(_mechtype, 51, 6);
  hoc_register_dparam_semantics(_mechtype, 0, "area");
  hoc_register_dparam_semantics(_mechtype, 1, "pntproc");
  hoc_register_dparam_semantics(_mechtype, 2, "netsend");
  hoc_register_dparam_semantics(_mechtype, 3, "watch");
  hoc_register_dparam_semantics(_mechtype, 4, "watch");
  hoc_register_dparam_semantics(_mechtype, 5, "cvodeieq");
 	hoc_register_cvode(_mechtype, _ode_count, _ode_map, _ode_spec, _ode_matsol);
 	hoc_register_tolerance(_mechtype, _hoc_state_tol, &_atollist);
 add_nrn_has_net_event(_mechtype);
 pnt_receive[_mechtype] = _net_receive;
 pnt_receive_init[_mechtype] = _net_init;
 pnt_receive_size[_mechtype] = 3;
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 AdExp C:/Users/teppola/Documents/OCNC_Project/Heidis_modeling_project_OCNC/data-analysis/Lonardoni2017_data/AdExp.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
static int _reset;
static char *modelname = "";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
 
#define _deriv1_advance _thread[0]._i
#define _dith1 1
#define _recurse _thread[2]._i
#define _newtonspace1 _thread[3]._pvoid
extern void* nrn_cons_newtonspace(int);
 
static int _ode_spec1(_threadargsproto_);
/*static int _ode_matsol1(_threadargsproto_);*/
 static int _slist2[8];
  static int _slist1[8], _dlist1[8];
 static int states(_threadargsproto_);
 
/*CVODE*/
 static int _ode_spec1 (double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {int _reset = 0; {
   Dy1 = - y1 / ( tDECAYnmda ) ;
   Dy2 = - y2 / ( tRISEnmda ) ;
   gNMDA = ( y2 - y1 ) ;
   iNMDA = ( gNMDA ) * ( vv - Erev ) * 1.0 / ( 1.0 + exp ( - ( vv - v0_block ) / k_block ) ) ;
   DgRISEexc = - gRISEexc / tRISE ;
   DgRISEinh = - gRISEinh / tRISE ;
   DgEXC = - gEXC / tauEXC ;
   DgINH = - gINH / tauINH ;
   gAMPA = ( gRISEexc - gEXC ) ;
   iAMPA = gAMPA * ( vv - erev ) ;
   gGABA = ( gRISEinh - gINH ) ;
   iGABA = - gGABA * ( vv - irev ) ;
   iTotal = - ww + iEXT + iAMPA + iGABA + iNMDA ;
   Dvv = 1.0 / C * ( - G_l * ( vv - E_l ) + G_l * Delta_T * exp ( ( ( vv - V_thre ) ) / Delta_T ) + iTotal ) ;
   Dww = 1.0 / tau_w * ( a * ( vv - E_l ) - ww ) ;
   }
 return _reset;
}
 static int _ode_matsol1 (double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {
 Dy1 = Dy1  / (1. - dt*( ( - 1.0 ) / ( tDECAYnmda ) )) ;
 Dy2 = Dy2  / (1. - dt*( ( - 1.0 ) / ( tRISEnmda ) )) ;
 gNMDA = ( y2 - y1 ) ;
 iNMDA = ( gNMDA ) * ( vv - Erev ) * 1.0 / ( 1.0 + exp ( - ( vv - v0_block ) / k_block ) ) ;
 DgRISEexc = DgRISEexc  / (1. - dt*( ( - 1.0 ) / tRISE )) ;
 DgRISEinh = DgRISEinh  / (1. - dt*( ( - 1.0 ) / tRISE )) ;
 DgEXC = DgEXC  / (1. - dt*( ( - 1.0 ) / tauEXC )) ;
 DgINH = DgINH  / (1. - dt*( ( - 1.0 ) / tauINH )) ;
 gAMPA = ( gRISEexc - gEXC ) ;
 iAMPA = gAMPA * ( vv - erev ) ;
 gGABA = ( gRISEinh - gINH ) ;
 iGABA = - gGABA * ( vv - irev ) ;
 iTotal = - ww + iEXT + iAMPA + iGABA + iNMDA ;
 Dvv = Dvv  / (1. - dt*( (( 1.0 / C * ( - G_l * ( ( vv  + .001) - E_l ) + G_l * Delta_T * exp ( ( ( ( vv  + .001) - V_thre ) ) / Delta_T ) + iTotal ) ) - ( 1.0 / C * ( - G_l * ( vv - E_l ) + G_l * Delta_T * exp ( ( ( vv - V_thre ) ) / Delta_T ) + iTotal )  )) / .001 )) ;
 Dww = Dww  / (1. - dt*( ( 1.0 / tau_w )*( ( ( - 1.0 ) ) ) )) ;
  return 0;
}
 /*END CVODE*/
 
static int states (double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {int _reset=0; int error = 0;
 { double* _savstate1 = _thread[_dith1]._pval;
 double* _dlist2 = _thread[_dith1]._pval + 8;
 int _counte = -1;
 if (!_recurse) {
 _recurse = 1;
 {int _id; for(_id=0; _id < 8; _id++) { _savstate1[_id] = _p[_slist1[_id]];}}
 error = nrn_newton_thread(_newtonspace1, 8,_slist2, _p, states, _dlist2, _ppvar, _thread, _nt);
 _recurse = 0; if(error) {abort_run(error);}}
 {
   Dy1 = - y1 / ( tDECAYnmda ) ;
   Dy2 = - y2 / ( tRISEnmda ) ;
   gNMDA = ( y2 - y1 ) ;
   iNMDA = ( gNMDA ) * ( vv - Erev ) * 1.0 / ( 1.0 + exp ( - ( vv - v0_block ) / k_block ) ) ;
   DgRISEexc = - gRISEexc / tRISE ;
   DgRISEinh = - gRISEinh / tRISE ;
   DgEXC = - gEXC / tauEXC ;
   DgINH = - gINH / tauINH ;
   gAMPA = ( gRISEexc - gEXC ) ;
   iAMPA = gAMPA * ( vv - erev ) ;
   gGABA = ( gRISEinh - gINH ) ;
   iGABA = - gGABA * ( vv - irev ) ;
   iTotal = - ww + iEXT + iAMPA + iGABA + iNMDA ;
   Dvv = 1.0 / C * ( - G_l * ( vv - E_l ) + G_l * Delta_T * exp ( ( ( vv - V_thre ) ) / Delta_T ) + iTotal ) ;
   Dww = 1.0 / tau_w * ( a * ( vv - E_l ) - ww ) ;
   {int _id; for(_id=0; _id < 8; _id++) {
if (_deriv1_advance) {
 _dlist2[++_counte] = _p[_dlist1[_id]] - (_p[_slist1[_id]] - _savstate1[_id])/dt;
 }else{
_dlist2[++_counte] = _p[_slist1[_id]] - _savstate1[_id];}}}
 } }
 return _reset;}
 
static double _watch1_cond(_pnt) Point_process* _pnt; {
 	double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
	_thread= (Datum*)0; _nt = (_NrnThread*)_pnt->_vnt;
 	_p = _pnt->_prop->param; _ppvar = _pnt->_prop->dparam;
	v = NODEV(_pnt->node);
	return  ( vv ) - ( V_spike ) ;
}
 
static void _net_receive (_pnt, _args, _lflag) Point_process* _pnt; double* _args; double _lflag; 
{  double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   int _watch_rm = 0;
   _thread = (Datum*)0; _nt = (_NrnThread*)_pnt->_vnt;   _p = _pnt->_prop->param; _ppvar = _pnt->_prop->dparam;
  if (_tsav > t){ extern char* hoc_object_name(); hoc_execerror(hoc_object_name(_pnt->ob), ":Event arrived out of order. Must call ParallelContext.set_maxstep AFTER assigning minimum NetCon.delay");}
 _tsav = t;   if (_lflag == 1. ) {*(_tqitem) = 0;}
 {
   if ( _lflag  == 1.0 ) {
       _nrn_watch_activate(_watch_array, _watch1_cond, 1, _pnt, _watch_rm++, 2.0);
 }
   else if ( _lflag  == 2.0 ) {
     net_event ( _pnt, t ) ;
       if (nrn_netrec_state_adjust && !cvode_active_){
    /* discon state adjustment for general derivimplicit and KINETIC case */
    int __i, __neq = 8;
    double __state = vv;
    double __primary_delta = (V_reset) - __state;
    double __dtsav = dt;
    for (__i = 0; __i < __neq; ++__i) {
      _p[_dlist1[__i]] = 0.0;
    }
    _p[_dlist1[6]] = __primary_delta;
    dt *= 0.5;
    v = NODEV(_pnt->node);
#if NRN_VECTORIZED
    _thread = _nt->_ml_list[_mechtype]->_thread;
#endif
    _ode_matsol_instance1(_threadargs_);
    dt = __dtsav;
    for (__i = 0; __i < __neq; ++__i) {
      _p[_slist1[__i]] += _p[_dlist1[__i]];
    }
  } else {
 vv = V_reset ;
       }
   if (nrn_netrec_state_adjust && !cvode_active_){
    /* discon state adjustment for general derivimplicit and KINETIC case */
    int __i, __neq = 8;
    double __state = ww;
    double __primary_delta = (ww + b) - __state;
    double __dtsav = dt;
    for (__i = 0; __i < __neq; ++__i) {
      _p[_dlist1[__i]] = 0.0;
    }
    _p[_dlist1[7]] = __primary_delta;
    dt *= 0.5;
    v = NODEV(_pnt->node);
#if NRN_VECTORIZED
    _thread = _nt->_ml_list[_mechtype]->_thread;
#endif
    _ode_matsol_instance1(_threadargs_);
    dt = __dtsav;
    for (__i = 0; __i < __neq; ++__i) {
      _p[_slist1[__i]] += _p[_dlist1[__i]];
    }
  } else {
 ww = ww + b ;
       }
 }
   else if ( _lflag  == 0.0 ) {
     if ( _args[0] > 0.0 ) {
       if ( _args[0] > 1000.0 ) {
           if (nrn_netrec_state_adjust && !cvode_active_){
    /* discon state adjustment for general derivimplicit and KINETIC case */
    int __i, __neq = 8;
    double __state = gEXC;
    double __primary_delta = (gEXC + rpeso) - __state;
    double __dtsav = dt;
    for (__i = 0; __i < __neq; ++__i) {
      _p[_dlist1[__i]] = 0.0;
    }
    _p[_dlist1[4]] = __primary_delta;
    dt *= 0.5;
    v = NODEV(_pnt->node);
#if NRN_VECTORIZED
    _thread = _nt->_ml_list[_mechtype]->_thread;
#endif
    _ode_matsol_instance1(_threadargs_);
    dt = __dtsav;
    for (__i = 0; __i < __neq; ++__i) {
      _p[_slist1[__i]] += _p[_dlist1[__i]];
    }
  } else {
 gEXC = gEXC + rpeso ;
           }
   if (nrn_netrec_state_adjust && !cvode_active_){
    /* discon state adjustment for general derivimplicit and KINETIC case */
    int __i, __neq = 8;
    double __state = gRISEexc;
    double __primary_delta = (gRISEexc + rpeso) - __state;
    double __dtsav = dt;
    for (__i = 0; __i < __neq; ++__i) {
      _p[_dlist1[__i]] = 0.0;
    }
    _p[_dlist1[2]] = __primary_delta;
    dt *= 0.5;
    v = NODEV(_pnt->node);
#if NRN_VECTORIZED
    _thread = _nt->_ml_list[_mechtype]->_thread;
#endif
    _ode_matsol_instance1(_threadargs_);
    dt = __dtsav;
    for (__i = 0; __i < __neq; ++__i) {
      _p[_slist1[__i]] += _p[_dlist1[__i]];
    }
  } else {
 gRISEexc = gRISEexc + rpeso ;
           }
 }
       else {
           if (nrn_netrec_state_adjust && !cvode_active_){
    /* discon state adjustment for general derivimplicit and KINETIC case */
    int __i, __neq = 8;
    double __state = y1;
    double __primary_delta = (y1 + ( _args[0] * ( mNMDA ) )) - __state;
    double __dtsav = dt;
    for (__i = 0; __i < __neq; ++__i) {
      _p[_dlist1[__i]] = 0.0;
    }
    _p[_dlist1[0]] = __primary_delta;
    dt *= 0.5;
    v = NODEV(_pnt->node);
#if NRN_VECTORIZED
    _thread = _nt->_ml_list[_mechtype]->_thread;
#endif
    _ode_matsol_instance1(_threadargs_);
    dt = __dtsav;
    for (__i = 0; __i < __neq; ++__i) {
      _p[_slist1[__i]] += _p[_dlist1[__i]];
    }
  } else {
 y1 = y1 + ( _args[0] * ( mNMDA ) ) ;
           }
   if (nrn_netrec_state_adjust && !cvode_active_){
    /* discon state adjustment for general derivimplicit and KINETIC case */
    int __i, __neq = 8;
    double __state = y2;
    double __primary_delta = (y2 + ( _args[0] * ( mNMDA ) )) - __state;
    double __dtsav = dt;
    for (__i = 0; __i < __neq; ++__i) {
      _p[_dlist1[__i]] = 0.0;
    }
    _p[_dlist1[1]] = __primary_delta;
    dt *= 0.5;
    v = NODEV(_pnt->node);
#if NRN_VECTORIZED
    _thread = _nt->_ml_list[_mechtype]->_thread;
#endif
    _ode_matsol_instance1(_threadargs_);
    dt = __dtsav;
    for (__i = 0; __i < __neq; ++__i) {
      _p[_slist1[__i]] += _p[_dlist1[__i]];
    }
  } else {
 y2 = y2 + ( _args[0] * ( mNMDA ) ) ;
           }
   if (nrn_netrec_state_adjust && !cvode_active_){
    /* discon state adjustment for general derivimplicit and KINETIC case */
    int __i, __neq = 8;
    double __state = gEXC;
    double __primary_delta = (gEXC + _args[0] * _args[1] * ( mAMPA )) - __state;
    double __dtsav = dt;
    for (__i = 0; __i < __neq; ++__i) {
      _p[_dlist1[__i]] = 0.0;
    }
    _p[_dlist1[4]] = __primary_delta;
    dt *= 0.5;
    v = NODEV(_pnt->node);
#if NRN_VECTORIZED
    _thread = _nt->_ml_list[_mechtype]->_thread;
#endif
    _ode_matsol_instance1(_threadargs_);
    dt = __dtsav;
    for (__i = 0; __i < __neq; ++__i) {
      _p[_slist1[__i]] += _p[_dlist1[__i]];
    }
  } else {
 gEXC = gEXC + _args[0] * _args[1] * ( mAMPA ) ;
           }
   if (nrn_netrec_state_adjust && !cvode_active_){
    /* discon state adjustment for general derivimplicit and KINETIC case */
    int __i, __neq = 8;
    double __state = gRISEexc;
    double __primary_delta = (gRISEexc + _args[0] * _args[1] * ( mAMPA )) - __state;
    double __dtsav = dt;
    for (__i = 0; __i < __neq; ++__i) {
      _p[_dlist1[__i]] = 0.0;
    }
    _p[_dlist1[2]] = __primary_delta;
    dt *= 0.5;
    v = NODEV(_pnt->node);
#if NRN_VECTORIZED
    _thread = _nt->_ml_list[_mechtype]->_thread;
#endif
    _ode_matsol_instance1(_threadargs_);
    dt = __dtsav;
    for (__i = 0; __i < __neq; ++__i) {
      _p[_slist1[__i]] += _p[_dlist1[__i]];
    }
  } else {
 gRISEexc = gRISEexc + _args[0] * _args[1] * ( mAMPA ) ;
           }
 }
       }
     else {
         if (nrn_netrec_state_adjust && !cvode_active_){
    /* discon state adjustment for general derivimplicit and KINETIC case */
    int __i, __neq = 8;
    double __state = gINH;
    double __primary_delta = (gINH + _args[0]) - __state;
    double __dtsav = dt;
    for (__i = 0; __i < __neq; ++__i) {
      _p[_dlist1[__i]] = 0.0;
    }
    _p[_dlist1[5]] = __primary_delta;
    dt *= 0.5;
    v = NODEV(_pnt->node);
#if NRN_VECTORIZED
    _thread = _nt->_ml_list[_mechtype]->_thread;
#endif
    _ode_matsol_instance1(_threadargs_);
    dt = __dtsav;
    for (__i = 0; __i < __neq; ++__i) {
      _p[_slist1[__i]] += _p[_dlist1[__i]];
    }
  } else {
 gINH = gINH + _args[0] ;
         }
   if (nrn_netrec_state_adjust && !cvode_active_){
    /* discon state adjustment for general derivimplicit and KINETIC case */
    int __i, __neq = 8;
    double __state = gRISEinh;
    double __primary_delta = (gRISEinh + _args[0]) - __state;
    double __dtsav = dt;
    for (__i = 0; __i < __neq; ++__i) {
      _p[_dlist1[__i]] = 0.0;
    }
    _p[_dlist1[3]] = __primary_delta;
    dt *= 0.5;
    v = NODEV(_pnt->node);
#if NRN_VECTORIZED
    _thread = _nt->_ml_list[_mechtype]->_thread;
#endif
    _ode_matsol_instance1(_threadargs_);
    dt = __dtsav;
    for (__i = 0; __i < __neq; ++__i) {
      _p[_slist1[__i]] += _p[_dlist1[__i]];
    }
  } else {
 gRISEinh = gRISEinh + _args[0] ;
         }
 }
     }
   } }
 
static void _net_init(Point_process* _pnt, double* _args, double _lflag) {
       double* _p = _pnt->_prop->param;
    Datum* _ppvar = _pnt->_prop->dparam;
    Datum* _thread = (Datum*)0;
    _NrnThread* _nt = (_NrnThread*)_pnt->_vnt;
 _args[1] = _args[1] ;
   _args[2] = _args[2] ;
   }
 
static int _ode_count(int _type){ return 8;}
 
static void _ode_spec(_NrnThread* _nt, _Memb_list* _ml, int _type) {
   double* _p; Datum* _ppvar; Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
     _ode_spec1 (_p, _ppvar, _thread, _nt);
 }}
 
static void _ode_map(int _ieq, double** _pv, double** _pvdot, double* _pp, Datum* _ppd, double* _atol, int _type) { 
	double* _p; Datum* _ppvar;
 	int _i; _p = _pp; _ppvar = _ppd;
	_cvode_ieq = _ieq;
	for (_i=0; _i < 8; ++_i) {
		_pv[_i] = _pp + _slist1[_i];  _pvdot[_i] = _pp + _dlist1[_i];
		_cvode_abstol(_atollist, _atol, _i);
	}
 }
 
static void _ode_matsol_instance1(_threadargsproto_) {
 _ode_matsol1 (_p, _ppvar, _thread, _nt);
 }
 
static void _ode_matsol(_NrnThread* _nt, _Memb_list* _ml, int _type) {
   double* _p; Datum* _ppvar; Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
 _ode_matsol_instance1(_threadargs_);
 }}
 
static void _thread_mem_init(Datum* _thread) {
   _thread[_dith1]._pval = (double*)ecalloc(16, sizeof(double));
   _newtonspace1 = nrn_cons_newtonspace(8);
 }
 
static void _thread_cleanup(Datum* _thread) {
   free((void*)(_thread[_dith1]._pval));
   nrn_destroy_newtonspace(_newtonspace1);
 }

static void initmodel(double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {
  int _i; double _save;{
  gINH = gINH0;
  gEXC = gEXC0;
  gRISEexc = gRISEexc0;
  gRISEinh = gRISEinh0;
  vv = vv0;
  ww = ww0;
  y2 = y20;
  y1 = y10;
 {
   gEXC = 0.0 ;
   gINH = 0.0 ;
   gRISEexc = 0.0 ;
   gRISEinh = 0.0 ;
   gAMPA = 0.0 ;
   gGABA = 0.0 ;
   ww = 0.0 ;
   vv = - 60.0 ;
   net_send ( _tqitem, (double*)0, _ppvar[1]._pvoid, t +  0.0 , 1.0 ) ;
   }
 
}
}

static void nrn_init(_NrnThread* _nt, _Memb_list* _ml, int _type){
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; double _v; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
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
 initmodel(_p, _ppvar, _thread, _nt);
}
}

static double _nrn_current(double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt, double _v){double _current=0.;v=_v;{
} return _current;
}

static void nrn_cur(_NrnThread* _nt, _Memb_list* _ml, int _type) {
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; int* _ni; double _rhs, _v; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 
}
 
}

static void nrn_jacob(_NrnThread* _nt, _Memb_list* _ml, int _type) {
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml];
#if CACHEVEC
  if (use_cachevec) {
	VEC_D(_ni[_iml]) += _g;
  }else
#endif
  {
     _nd = _ml->_nodelist[_iml];
	NODED(_nd) += _g;
  }
 
}
 
}

static void nrn_state(_NrnThread* _nt, _Memb_list* _ml, int _type) {
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; double _v = 0.0; int* _ni; int _iml, _cntml;
double _dtsav = dt;
if (secondorder) { dt *= 0.5; }
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
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
 {  _deriv1_advance = 1;
 derivimplicit_thread(8, _slist1, _dlist1, _p, states, _ppvar, _thread, _nt);
_deriv1_advance = 0;
     if (secondorder) {
    int _i;
    for (_i = 0; _i < 8; ++_i) {
      _p[_slist1[_i]] += dt*_p[_dlist1[_i]];
    }}
 } {
   }
}}
 dt = _dtsav;
}

static void terminal(){}

static void _initlists(){
 double _x; double* _p = &_x;
 int _i; static int _first = 1;
  if (!_first) return;
 _slist1[0] = &(y1) - _p;  _dlist1[0] = &(Dy1) - _p;
 _slist1[1] = &(y2) - _p;  _dlist1[1] = &(Dy2) - _p;
 _slist1[2] = &(gRISEexc) - _p;  _dlist1[2] = &(DgRISEexc) - _p;
 _slist1[3] = &(gRISEinh) - _p;  _dlist1[3] = &(DgRISEinh) - _p;
 _slist1[4] = &(gEXC) - _p;  _dlist1[4] = &(DgEXC) - _p;
 _slist1[5] = &(gINH) - _p;  _dlist1[5] = &(DgINH) - _p;
 _slist1[6] = &(vv) - _p;  _dlist1[6] = &(Dvv) - _p;
 _slist1[7] = &(ww) - _p;  _dlist1[7] = &(Dww) - _p;
 _slist2[0] = &(gINH) - _p;
 _slist2[1] = &(gEXC) - _p;
 _slist2[2] = &(gRISEexc) - _p;
 _slist2[3] = &(gRISEinh) - _p;
 _slist2[4] = &(vv) - _p;
 _slist2[5] = &(ww) - _p;
 _slist2[6] = &(y2) - _p;
 _slist2[7] = &(y1) - _p;
_first = 0;
}

#if defined(__cplusplus)
} /* extern "C" */
#endif

#if NMODL_TEXT
static const char* nmodl_filename = "AdExp.mod";
static const char* nmodl_file_text = 
  ": AdExp IF (Adaptive Exponential Intgrate and Fire)\n"
  ": All printfs are for debugging purposes and can be ignored\n"
  "\n"
  "\n"
  "NEURON {\n"
  "\n"
  "    POINT_PROCESS AdExp\n"
  "	RANGE num\n"
  "    RANGE G_l, Delta_T, tau_w, a, b, E_l, I_ext, C\n"
  "    RANGE taue, taui, erev, irev\n"
  "    RANGE V_reset, V_thre, tRISE,tot\n"
  "    RANGE gRISEinh,gRISEexc, iEXT,tauINH,tauEXC,tRISE,tRISEnmda,tDECAYnmda,gAMPA,gGABA\n"
  "    RANGE iNMDA,iAMPA,iGABA,iTotal,gGABA,gAMPA,gNMDA,mAMPA,mNMDA,rpeso,Erev,v0_block,k_block,tau_r, tau_d\n"
  "\n"
  "}\n"
  "\n"
  "UNITS {\n"
  "\n"
  "	(mV) = (millivolt)\n"
  "	(pA) = (picoamp)\n"
  "	(uS) = (microsiemens)\n"
  "    (nS) = (nanosiemens)\n"
  "	(pS) = (picosiemens)\n"
  "}\n"
  "\n"
  "PARAMETER {\n"
  "	num = -1  					:Cell ID\n"
  "	V_reset = -60	(mV)		:reset potential after spike\n"
  "	V_thre = -50.4	(mV)		:threshold for spike detection\n"
  "	V_spike = -10	(mV)	    :value of spike\n"
  "	a = 4			(nS)		:coupling with adaptive variable\n"
  "	b = 80.5		(pA)		:adaptive increment\n"
  "	tau_w = 144		(ms)		:adaptive time costant\n"
  "	E_l = -70.6		(mV)		:resting potential for leak term\n"
  "	G_l = 30		(nS)	    :leak conduptance\n"
  "	Delta_T = 2		(mV)	    :speed of exp\n"
  "\n"
  "	eps = 0.1		(ms) 		:smaller than time step\n"
  "	iEXT = 0        (pA) 	    :external current\n"
  "	C = 281 		(nS)		:membrane conduptance\n"
  "    tauEXC = 3 		(ms) 		:excitatory time costant\n"
  "    tauINH = 8		(ms) 		:inhibitory time costant\n"
  "    erev = 0        (mV) 		:reverse potential for exc\n"
  "    irev = -70		(mV)		:reverse potential for inh\n"
  "    tRISE = 1 		(ms) 		:time to activate the synapse\n"
  "\n"
  "    mNMDA=.05\n"
  "    mAMPA=1\n"
  "    v0_block 	= -50\n"
  "	k_block 	= 8\n"
  "\n"
  "	tRISEnmda = 5.63    (ms)        :Chapman DE 2003, Table 1 - rise tau\n"
  "	tDECAYnmda = 140     (ms)        :Chapman DE 2003, Fig 2B Ventromedial - decay tau\n"
  "\n"
  "	Erev= 0         (mV)        :reversal potential, Dalby 2003\n"
  "	rpeso=12        (nS)\n"
  "\n"
  "}\n"
  "\n"
  "ASSIGNED {\n"
  "    gNMDA           (nS)\n"
  "    gAMPA           (nS)\n"
  "    gGABA           (nS)\n"
  "	iAMPA           (pA)\n"
  "    iGABA           (pA)\n"
  "	iNMDA           (pA)\n"
  "    iNOISE          (pA)\n"
  "    iTotal          (pA)\n"
  "\n"
  "\n"
  "}\n"
  "\n"
  "STATE {\n"
  "    vv              (mV)\n"
  "    ww              (pA)\n"
  "    gEXC            (nS)\n"
  "    gINH            (nS)\n"
  "    gRISEexc        (nS)\n"
  "    gRISEinh        (nS)\n"
  "    y1              (nS)\n"
  "    y2              (nS)\n"
  "}\n"
  "\n"
  "INITIAL { :initializion and activation of the process\n"
  "\n"
  "    gEXC=0          (nS)\n"
  "    gINH=0          (nS)\n"
  "    gRISEexc=0      (nS)\n"
  "    gRISEinh=0      (nS)\n"
  "    gAMPA=0         (nS)\n"
  "    gGABA=0         (nS)\n"
  "    ww=0            (pA)\n"
  "    vv=-60          (mV)\n"
  "    net_send(0,1)\n"
  "}\n"
  "\n"
  "\n"
  "BREAKPOINT {\n"
  "\n"
  "\n"
  "\n"
  "SOLVE states METHOD derivimplicit\n"
  "\n"
  "\n"
  "}\n"
  "\n"
  "DERIVATIVE states {\n"
  "\n"
  "	y1' = -y1 / (tDECAYnmda)\n"
  "	y2' = -y2 / (tRISEnmda)\n"
  "    gNMDA = (y2-y1)\n"
  "	iNMDA = (gNMDA)* (vv - Erev)*1 / ( 1 + exp ( - ( vv - v0_block ) / k_block ))\n"
  "\n"
  "    gRISEexc' = - gRISEexc/tRISE        :rise of excitatory synapse\n"
  "    gRISEinh'= - gRISEinh/tRISE         :rise of inhibitory synapse\n"
  "    gEXC' = - gEXC /tauEXC              :decay of exc conduptance\n"
  "    gINH' = - gINH /tauINH              :decay of inh conduptance\n"
  "\n"
  "    gAMPA = (gRISEexc-gEXC)\n"
  "	iAMPA = gAMPA * ( vv - erev )\n"
  "\n"
  "    gGABA = (gRISEinh-gINH)\n"
  "    iGABA = -gGABA * ( vv - irev )\n"
  "\n"
  "\n"
  "    														:iNOISE = - (x2-x1) * ( vv - erev )\n"
  "    iTotal=-ww+iEXT+iAMPA+iGABA+iNMDA  						:+iNOISE\n"
  "\n"
  "    vv'=1/C*(-G_l*(vv-E_l)+G_l*Delta_T*exp(((vv-V_thre))/Delta_T)+iTotal)\n"
  "    ww' = 1/tau_w*(a*(vv-E_l)-ww)\n"
  "\n"
  "    }\n"
  "\n"
  "\n"
  "NET_RECEIVE (w0,w1,w2) {\n"
  "    INITIAL {\n"
  "        w1=w1\n"
  "        w2=w2\n"
  "        }\n"
  "    if (flag == 1) {                            :wait for membrane potential to reach threshold\n"
  "        WATCH (vv>V_spike) 2\n"
  "    }\n"
  "    else if (flag == 2) {                       :threshold reached, fire and reset the model\n"
  "\n"
  "    net_event(t)\n"
  "    vv = V_reset\n"
  "    ww = ww+b\n"
  "\n"
  "    }\n"
  "    else if (flag==0) {                         :an event has arrived -> synaptic activation\n"
  "\n"
  "    if ( w0 > 0 ) {                             :positive weight -> excitatory synapse\n"
  "\n"
  "        if (w0>1000) {                              :overlimit conductance -> poisson process\n"
  "            gEXC=gEXC+rpeso\n"
  "            gRISEexc=gRISEexc+rpeso\n"
  "            :x1=x1+rpeso\n"
  "            :x2=x2+rpeso\n"
  "        }\n"
  "\n"
  "        else {                                      :genuine excitation -> update synconductance\n"
  "\n"
  "            y1=y1+(w0 *(mNMDA))\n"
  "            y2=y2+(w0 *(mNMDA))\n"
  "            gEXC=gEXC+w0 *w1*(mAMPA)\n"
  "            gRISEexc=gRISEexc+w0 *w1*(mAMPA)\n"
  "        }\n"
  "\n"
  "    }\n"
  "    else {                                          :genuine inhibition -> update synconductance\n"
  "\n"
  "        gINH=gINH+w0\n"
  "        gRISEinh=gRISEinh+w0\n"
  "\n"
  "        }\n"
  "    }\n"
  "}\n"
  ;
#endif
