// -*- C++ -*-
// 
//  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// 
//                               Charles A. Williams
//                        Rensselaer Polytechnic Institute
//                        (C) 2004 All Rights Reserved
// 
//  Copyright 2004 Rensselaer Polytechnic Institute.
//  All worldwide rights reserved.  A license to use, copy, modify and
//  distribute this software for non-commercial research purposes only
//  is hereby granted, provided that this copyright notice and
//  accompanying disclaimer is not modified or removed from the software.
//
//  DISCLAIMER:  The software is distributed "AS IS" without any express
//  or implied warranty, including but not limited to, any implied
//  warranties of merchantability or fitness for a particular purpose
//  or any warranty of non-infringement of any current or pending patent
//  rights.  The authors of the software make no representations about
//  the suitability of this software for any particular purpose.  The
//  entire risk as to the quality and performance of the software is with
//  the user.  Should the software prove defective, the user assumes the
//  cost of all necessary servicing, repair or correction.  In
//  particular, neither Rensselaer Polytechnic Institute, nor the authors
//  of the software are liable for any indirect, special, consequential,
//  or incidental damages related to the software, to the maximum extent
//  the law permits.
// 
//  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// 

#if !defined(lithomop3d_externs_h)
#define lithomop3d_externs_h

#if defined(NEEDS_F77_TRANSLATION)

#if defined(F77EXTERNS_LOWERCASE_TRAILINGBAR)

#define autoprestr_f autoprestr_
#define cmp_stiffsz_f cmp_stiffsz_
#define create_id_f create_id_
#define elastc_f elastc_
#define id_split_f id_split_
#define lnklst_f lnklst_
#define local_f local_
#define localf_f localf_
#define localx_f localx_
#define makemsr_f makemsr_
#define matmod_def_f matmod_def_
#define nfind_f nfind_
#define preshape_f preshape_
#define read_bc_f read_bc_
#define read_connect_f read_connect_
#define read_coords_f read_coords_
#define read_diff_f read_diff_
#define read_fuldat_f read_fuldat_
#define read_hist_f read_hist_
#define read_mathist_f read_mathist_
// #define read_prestr_f read_prestr_
#define read_skew_f read_skew_
#define read_slip_f read_slip_
#define read_split_f read_split_
#define read_stateout_f read_stateout_
#define read_timdat_f read_timdat_
// #define read_traction_f read_traction_
#define read_wink_f read_wink_
#define read_winkx_f read_winkx_
#define scan_bc_f scan_bc_
#define scan_connect_f scan_connect_
#define scan_coords_f scan_coords_
#define scan_diff_f scan_diff_
#define scan_fuldat_f scan_fuldat_
#define scan_hist_f scan_hist_
// #define scan_prestr_f scan_prestr_
#define scan_skew_f scan_skew_
#define scan_slip_f scan_slip_
#define scan_split_f scan_split_
#define scan_timdat_f scan_timdat_
// #define scan_traction_f scan_traction_
#define scan_wink_f scan_wink_
#define scan_winkx_f scan_winkx_
#define viscos_f viscos_
#define write_element_info_f write_element_info_
#define write_global_info_f write_global_info_
#define write_props_f write_props_
#define write_sparse_info_f write_sparse_info_
#define write_strscomp_f write_strscomp_
#define write_subiter_f write_subiter_
#define write_ucd_mesh_f write_ucd_mesh_

#elif defined(F77EXTERNS_NOTRAILINGBAR)

#define autoprestr_f autoprestr
#define cmp_stiffsz_f cmp_stiffsz
#define create_id_f create_id
#define elastc_f elastc
#define id_split_f id_split
#define lnklst_f lnklst
#define local_f local
#define localf_f localf
#define localx_f localx
#define makemsr_f makemsr
#define matmod_def_f matmod_def
#define nfind_f nfind
#define preshape_f preshape
#define read_bc_f read_bc
#define read_connect_f read_connect
#define read_coords_f read_coords
#define read_diff_f read_diff
#define read_fuldat_f read_fuldat
#define read_hist_f read_hist
#define read_mathist_f read_mathist
// #define read_prestr_f read_prestr
#define read_skew_f read_skew
#define read_slip_f read_slip
#define read_split_f read_split
#define read_stateout_f read_stateout
#define read_timdat_f read_timdat
// #define read_traction_f read_traction
#define read_wink_f read_wink
#define read_winkx_f read_winkx
#define scan_bc_f scan_bc
#define scan_connect_f scan_connect
#define scan_coords_f scan_coords
#define scan_diff_f scan_diff
#define scan_fuldat_f scan_fuldat
#define scan_hist_f scan_hist
// #define scan_prestr_f scan_prestr
#define scan_skew_f scan_skew
#define scan_slip_f scan_slip
#define scan_split_f scan_split
#define scan_timdat_f scan_timdat
// #define scan_traction_f scan_traction
#define scan_wink_f scan_wink
#define scan_winkx_f scan_winkx
#define viscos_f viscos
#define write_element_info_f write_element_info
#define write_global_info_f write_global_info
#define write_props_f write_props
#define write_sparse_info_f write_sparse_info
#define write_strscomp_f write_strscomp
#define write_subiter_f write_subiter
#define write_ucd_mesh_f write_ucd_mesh

#elif defined(F77EXTERNS_EXTRATRAILINGBAR)

#define autoprestr_f autoprestr__
#define cmp_stiffsz_f cmp_stiffsz__
#define create_id_f create_id__
#define elastc_f elastc__
#define id_split_f id_split__
#define lnklst_f lnklst__
#define local_f local__
#define localf_f localf__
#define localx_f localx__
#define makemsr_f makemsr__
#define matmod_def_f matmod_def__
#define nfind_f nfind__
#define preshape_f preshape__
#define read_bc_f read_bc__
#define read_connect_f read_connect__
#define read_coords_f read_coords__
#define read_diff_f read_diff__
#define read_fuldat_f read_fuldat__
#define read_hist_f read_hist__
#define read_mathist_f read_mathist__
// #define read_prestr_f read_prestr__
#define read_skew_f read_skew__
#define read_slip_f read_slip__
#define read_split_f read_split__
#define read_stateout_f read_stateout__
#define read_timdat_f read_timdat__
// #define read_traction_f read_traction__
#define read_wink_f read_wink__
#define read_winkx_f read_winkx__
#define scan_bc_f scan_bc__
#define scan_connect_f scan_connect__
#define scan_coords_f scan_coords__
#define scan_diff_f scan_diff__
#define scan_fuldat_f scan_fuldat__
#define scan_hist_f scan_hist__
// #define scan_prestr_f scan_prestr__
#define scan_skew_f scan_skew__
#define scan_slip_f scan_slip__
#define scan_split_f scan_split__
#define scan_timdat_f scan_timdat__
// #define scan_traction_f scan_traction__
#define scan_wink_f scan_wink__
#define scan_winkx_f scan_winkx__
#define viscos_f viscos__
#define write_element_info_f write_element_info__
#define write_global_info_f write_global_info__
#define write_props_f write_props__
#define write_sparse_info_f write_sparse_info__
#define write_strscomp_f write_strscomp__
#define write_subiter_f write_subiter__
#define write_ucd_mesh_f write_ucd_mesh__

#elif defined(F77EXTERNS_UPPERCASE_NOTRAILINGBAR)

#define autoprestr_f AUTOPRESTR
#define cmp_stiffsz_f CMP_STIFFSZ
#define create_id_f CREATE_ID
#define elastc_f ELASTC
#define id_split_f ID_SPLIT
#define lnklst_f LNKLST
#define local_f LOCAL
#define localf_f LOCALF
#define localx_f LOCALX
#define makemsr_f MAKEMSR
#define matmod_def_f MATMOD_DEF
#define nfind_f NFIND
#define preshape_f PRESHAPE
#define read_bc_f READ_BC
#define read_connect_f READ_CONNECT
#define read_coords_f READ_COORDS
#define read_diff_f READ_DIFF
#define read_fuldat_f READ_FULDAT
#define read_hist_f READ_HIST
#define read_mathist_f READ_MATHIST
// #define read_prestr_f READ_PRESTR
#define read_skew_f READ_SKEW
#define read_slip_f READ_SLIP
#define read_split_f READ_SPLIT
#define read_stateout_f READ_STATEOUT
#define read_timdat_f READ_TIMDAT
// #define read_traction_f READ_TRACTION
#define read_wink_f READ_WINK
#define read_winkx_f READ_WINKX
#define scan_bc_f SCAN_BC
#define scan_connect_f SCAN_CONNECT
#define scan_coords_f SCAN_COORDS
#define scan_diff_f SCAN_DIFF
#define scan_fuldat_f SCAN_FULDAT
#define scan_hist_f SCAN_HIST
// #define scan_prestr_f SCAN_PRESTR
#define scan_skew_f SCAN_SKEW
#define scan_slip_f SCAN_SLIP
#define scan_split_f SCAN_SPLIT
#define scan_timdat_f SCAN_TIMDAT
// #define scan_traction_f SCAN_TRACTION
#define scan_wink_f SCAN_WINK
#define scan_winkx_f SCAN_WINKX
#define viscos_f VISCOS
#define write_element_info_f WRITE_ELEMENT_INFO
#define write_global_info_f WRITE_GLOBAL_INFO
#define write_props_f WRITE_PROPS
#define write_sparse_info_f WRITE_SPARSE_INFO
#define write_strscomp_f WRITE_STRSCOMP
#define write_subiter_f WRITE_SUBITER
#define write_ucd_mesh_f WRITE_UCD_MESH

#elif defined(F77EXTERNS_COMPAQ_F90) || defined (F77EXTERNS_SOMETIMES_EXTRATRAILINGBAR)

// symbols that contain underbars get two underbars at the end
// symbols that do not contain underbars get one underbar at the end
// this applies to the FORTRAN external, not the local macro alias!!!

#define autoprestr_f autoprestr_
#define cmp_stiffsz_f cmp_stiffsz__
#define create_id_f create_id__
#define elastc_f elastc_
#define id_split_f id_split__
#define lnklst_f lnklst_
#define local_f local_
#define localf_f localf_
#define localx_f localx_
#define makemsr_f makemsr_
#define matmod_def_f matmod_def__
#define nfind_f nfind_
#define preshape_f preshape_
#define read_bc_f read_bc__
#define read_connect_f read_connect__
#define read_coords_f read_coords__
#define read_diff_f read_diff__
#define read_fuldat_f read_fuldat__
#define read_hist_f read_hist__
#define read_mathist_f read_mathist__
// #define read_prestr_f read_prestr__
#define read_skew_f read_skew__
#define read_slip_f read_slip__
#define read_split_f read_split__
#define read_stateout_f read_stateout__
#define read_timdat_f read_timdat__
// #define read_traction_f read_traction__
#define read_wink_f read_wink__
#define read_winkx_f read_winkx__
#define scan_bc_f scan_bc__
#define scan_connect_f scan_connect__
#define scan_coords_f scan_coords__
#define scan_diff_f scan_diff__
#define scan_fuldat_f scan_fuldat__
#define scan_hist_f scan_hist__
// #define scan_prestr_f scan_prestr__
#define scan_skew_f scan_skew__
#define scan_slip_f scan_slip__
#define scan_split_f scan_split__
#define scan_timdat_f scan_timdat__
// #define scan_traction_f scan_traction__
#define scan_wink_f scan_wink__
#define scan_winkx_f scan_winkx__
#define viscos_f viscos_
#define write_element_info_f write_element_info__
#define write_global_info_f write_global_info__
#define write_props_f write_props__
#define write_sparse_info_f write_sparse_info__
#define write_strscomp_f write_strscomp__
#define write_subiter_f write_subiter__
#define write_ucd_mesh_f write_ucd_mesh__

#else
#error Unknown translation for FORTRAN external symbols
#endif

#endif

#include <petscmat.h>

extern "C" {

  void autoprestr_f(Mat *,             // sparse
		    double *,
		    double *,
		    double *,
		    double *,          // force
		    double *,
		    double *,
		    double *,
		    double *,
		    double *,
		    double *,
		    double *,
		    int *,
		    double *,
		    double *,          // global
		    double *,
		    double *,
		    double *,
		    int *,
		    int *,
		    double *,
		    int *,
		    int *,             // bc
		    double *,
		    double *,          // slip
		    double *,
		    double *,
		    double *,
		    int *,
		    int *,
		    double *,
		    int *,
		    int *,
		    int *,
		    double *,          // split
		    int *,
		    double *,
		    double *,
		    double *,          // stiff
		    double *,
		    double *,          // element
		    double *,
		    double *,
		    double *,
		    int *,
		    int *,
		    int *,
		    int *,
		    int *,
		    int *,
		    int *,
		    int *,             // traction
		    int *,
		    int *,
		    double *,
		    double *,
		    double *,          // material
		    int *,
		    int *,
		    int *,
		    int *,
		    double *,          // element type
		    double *,
		    double *,
		    int *,
		    double *,          // time data
		    double *,
		    int *,
		    int *,
		    int *,
		    double *,
		    double *,
		    int *,
		    int *,
		    int *,
		    double *,
		    double *,
		    double *,
		    int *,
		    double *,          // iterations
		    double *,
		    double *,
		    double *,
		    double *,
		    double *,
		    double *,
		    int *,
		    double *,          // skew
		    int *,             // i/o info
		    int *,
		    int *,
		    int *,
		    char *,            // files
		    char *,
		    char *,
                    int *,
		    int *,             // error codes
		    char *,
		    int,               // string lengths
		    int,
		    int,
		    int);

  void cmp_stiffsz_f(int *,
	             int *,
	             int *,
	             int *,
	             int *,
	             int *,
	             int *,
	             int *,
	             int *,
	             int *,
	             char *,int);

  void create_id_f(int *,
		   int *,
		   int *,
		   int *,
		   int *,
		   int *,
		   int *,
		   int *,
		   int *);

  void elastc_f(Mat *,             // sparse
		double *,
		double *,
		double *,
		double *,          // force
		double *,
		double *,
		double *,
		double *,
		double *,
		double *,
		double *,
		int *,
		double *,
		double *,          // global
		double *,
		double *,
		double *,
		int *,
		int *,
		double *,
		int *,
		int *,             // bc
		double *,
		double *,          // slip
		double *,
		double *,
		double *,
		int *,
		int *,
		double *,
		int *,
		int *,
		int *,
		double *,          // split
		int *,
		double *,
		double *,
		double *,          // stiff
		double *,
		double *,          // element
		double *,
		double *,
		double *,
		int *,
		int *,
		int *,
		int *,
		int *,
		int *,
		int *,
		int *,             // traction
		int *,
		int *,
		double *,
		double *,
		double *,          // material
		int *,
		int *,
		int *,
		int *,
		double *,          // element type
		double *,
		double *,
		int *,
		double *,          // time data
		double *,
		int *,
		int *,
		int *,
		double *,
		double *,
		int *,
		int *,
		int *,
		double *,
		double *,
		double *,
		int *,
		double *,          // iterations
		double *,
		double *,
		double *,
		double *,
		double *,
		double *,
		int *,
		double *,          // skew
		int *,             // i/o info
		int *,
		int *,
		int *,
		char *,            // files
		char *,
		char *,
                int *,
                int *,
		int *,             // error codes
		char *,
		int,               // string lengths
		int,
		int,
		int);

  void id_split_f(int *,
		  int *,
		  int *,
		  int *,
		  int *,
		  int *,
		  int *,
		  char *,int);

  void lnklst_f(int *,
		int *,
		int *,
		int *,
		int *,
		int *,
		int *,
		int *,
		int *,
		int *,
		int *,
		int *,
		int *,
		int *,
		int *,
		char *,
		int);

  void local_f(int *,
	       int *,
	       int *,
	       int *,
	       int *,
	       int *,
	       int *,
	       int *);

  void localf_f(int *,
		int *,
		int *,
		int *,
		int *,
		int *,
		int *,
		int *);

  void localx_f(int *,
		int *,
		int *,
		int *,
		int *,
		int *,
		int *,
		int *,
		int *,
		int *);

  void makemsr_f(Mat *,
		 int *,
		 int *,
		 int *,
		 int *,
		 int *,
		 int *,
		 int *,
		 int *,
		 double *);

  void matmod_def_f(int *,
		    int *);

  void nfind_f(double *,
	       double *,
	       int *,
	       int *,
	       int *,
	       int *,
	       int *,
	       int *,
	       int *,
	       int *,
	       int *);

  void preshape_f(double *,
		  double *,
		  double *,
		  int *,
		  int *);

  void read_bc_f(double*,
		 double*,
		 double*,
		 double*,
		 int *,
		 int *,
		 int *,
		 int *,
		 int *,
		 int *,
		 int *,
		 char *,
		 char *,
		 int *,
		 char *,
		 int,
		 int,
		 int);

  void read_connect_f(int *,
		      int *,
		      int *,
		      int *,
		      int *,
		      int *,
		      int *,
		      int *,
		      int *,
		      int *,
		      int *,
		      int *,
		      int *,
		      int *,
		      int *,
		      int *,
		      int *,
		      int *,
		      int *,
		      int *,
		      int *,
		      int *,
		      char *,
		      char *,
		      char *,
		      int *,
		      char *,
		      int,
		      int,
		      int,
		      int);

  void read_coords_f(double*,
		     double*,
		     int *,
		     int *,
		     int *,
		     int *,
		     int *,
		     int *,
		     char *,
		     char *,
		     char *,
		     int *,
		     char *,
		     int,
		     int,
		     int,
		     int);

  void read_diff_f(double*,
		   int *,
		   int *,
		   int *,
		   int *,
		   int *,
		   int *,
		   int *,
		   int *,
		   char *,
		   char *,
		   int *,
		   char *,
		   int,
		   int,
		   int);

  void read_fuldat_f(int *,
		     int *,
		     int *,
		     int *,
		     int *,
		     int *,
		     int *,
		     int *,
		     int *,
		     int *,
		     char *,
		     char *,
		     char *,
		     int *,
		     char *,
		     int,
		     int,
		     int,
		     int);

  void read_hist_f(double *,
		   double *,
		   int *,
		   int *,
		   int *,
		   int *,
		   int *,
		   char *,
		   char *,
		   int *,
		   char *,
		   int,
		   int,
		   int);

  void read_mathist_f(int *,
		      int *,
		      int *,
		      int *,
		      int *,
		      int *,
		      int *,
		      int *,
		      int *,
		      int *,
		      char *,
		      char *,
		      char *,
		      int *,
		      char *,
		      int,
		      int,
		      int,
		      int);

  // void read_prestr_f(double *,
		     // double *,
		     // double *,
		     // int *,
		     // int *,
		     // int *,
		     // int *,
		     // int *,
		     // int *,
		     // int *,
		     // int *,
		     // int *,
		     // int *,
		     // char *,
		     // char *,
		     // int,
		     // int);

  void read_skew_f(double *,
		   double *,
		   int *,
		   int *,
		   int *,
		   int *,
		   int *,
		   int *,
		   char *,
		   char *,
		   int *,
		   char *,
		   int,
		   int,
		   int);

  void read_slip_f(int *,
		   int *,
		   int *,
		   int *,
		   int *,
		   int *,
		   int *,
		   int *,
		   int *,
		   int *,
		   char *,
		   char *,
		   char *,
		   int *,
		   char *,
		   int,
		   int,
		   int,
		   int);

  void read_split_f(double *,
		    int *,
		    int *,
		    int *,
		    int *,
		    int *,
		    int *,
		    int *,
		    int *,
		    int *,
		    char *,
		    char *,
		    char *,
		    int *,
		    char *,
		    int,
		    int,
		    int,
		    int);

  void read_stateout_f(int *,
		       int *,
		       int *,
		       int *,
		       int *,
		       int *,
		       char *,
		       char *,
		       char *,
		       int *,
		       char *,
		       int,
		       int,
		       int,
		       int);

  void read_timdat_f(double *,
		     double *,
		     double *,
		     double *,
		     double *,
		     double *,
		     double *,
		     int *,
		     int *,
		     int *,
		     int *,
		     int *,
		     int *,
		     int *,
		     int *,
		     int *,
		     int *,
		     char *,
		     char *,
		     int *,
		     char *,
		     int,
		     int,
		     int);

  // void read_traction_f(double *,
		       // double *,
		       // double *,
		       // int *,
		       // int *,
		       // int *,
		       // int *,
		       // int *,
		       // int *,
		       // int *,
		       // int *,
		       // int *,
		       // int *,
		       // char *,
		       // char *,
		       // int,
		       // int);

  void read_wink_f(double *,
		   double *,
		   int *,
		   int *,
		   int *,
		   int *,
		   int *,
		   int *,
		   int *,
		   int *,
		   char *,
		   char *,
		   int *,
		   char *,
		   int,
		   int,
		   int);

  void read_winkx_f(double *,
		    double *,
		    int *,
		    int *,
		    int *,
		    int *,
		    int *,
		    int *,
		    int *,
		    int *,
		    char *,
		    char *,
		    int *,
		    char *,
		    int,
		    int,
		    int);

  void scan_bc_f(int *,
		 int *,
		 char *,
		 char *,
		 char *,
		 char *,
		 int *,
		 char *,
		 int,
		 int,
		 int,
		 int,
		 int);

  void scan_connect_f(int *,
		      int *,
		      int *,
		      int *,
		      int *,
		      int *,
		      int *,
		      int *,
		      int *,
		      char *,
		      int *,
		      char *,
		      int,
		      int);

  void scan_coords_f(int *,
		     int *,
		     char *,
		     char *,
		     int *,
		     char *,
		     int,
		     int,
		     int);

  void scan_diff_f(int *,
		   int *,
		   int *,
		   char *,
		   int *,
		   char *,
		   int,
		   int);

  void scan_fuldat_f(int *,
		     int *,
		     int *,
		     int *,
		     char *,
		     int *,
		     char *,
		     int,
		     int);

  void scan_hist_f(int *,
		   int *,
		   char *,
		   int *,
		   char *,
		   int,
		   int);

  // void scan_prestr_f(int *,
		     // int *,
		     // int *,
		     // int *,
		     // int *,
		     // int *,
		     // int *,
		     // char *,int);

  void scan_skew_f(int *,
		   int *,
		   char *,
		   char *,
		   int *,
		   char *,
		   int,
		   int,
		   int);

  void scan_slip_f(int *,
		   int *,
		   char *,
		   int *,
		   char *,
		   int,
		   int);

  void scan_split_f(int *,
		    int *,
		    char *,
		    int *,
		    char *,
		    int,
		    int);

  void scan_timdat_f(int *,
		     int *,
		     int *,
		     char *,
		     char *,
		     int *,
		     char *,
		     int,
		     int,
		     int);

  // void scan_traction_f(int *,
		       // int *,
		       // int *,
		       // int *,
		       // int *,
		       // char *,
		       // char *,
		       // int,
		       // int);

  void scan_wink_f(int *,
		   int *,
		   int *,
		   char *,
		   int *,
		   char *,
		   int,
		   int);

  void scan_winkx_f(int *,
		    int *,
		    int *,
		    int *,
		    char *,
		    int *,
		    char *,
		    int,
		    int);

  void viscos_f(Mat *,             // sparse
		double *,
		double *,
		double *,
		double *,          // force
		double *,
		double *,
		double *,
		double *,
		double *,
		double *,
		double *,
		int *,
		double *,
		double *,          // global
		double *,
		double *,
		double *,
		int *,
		int *,
		double *,
		int *,
		int *,             // bc
		double *,
		double *,          // slip
		double *,
		double *,
		double *,
		int *,
		int *,
		double *,
		int *,
		int *,
		int *,
		double *,          // split
		int *,
		double *,
		double *,
		double *,          // stiff
		double *,
		double *,          // element
		double *,
		double *,
		double *,
		int *,
		int *,
		int *,
		int *,
		int *,
		int *,
		int *,
		int *,             // traction
		int *,
		int *,
		double *,
		double *,
		double *,          // material
		int *,
		int *,
		int *,
		int *,
		double *,          // element type
		double *,
		double *,
		int *,
		double *,          // time data
		double *,
		int *,
		int *,
		int *,
		double *,
		double *,
		int *,
		int *,
		int *,
		double *,
		double *,
		double *,
		int *,
		double *,          // iterations
		double *,
		double *,
		double *,
		double *,
		double *,
		double *,
		int *,
		double *,          // skew
		int *,             // i/o info
		int *,
		int *,
		int *,
		int *,
		char *,            // files
		char *,
		char *,
                int *,
                int *,
		int *,             // error codes
		char *,
		int,               // string lengths
		int,
		int,
		int);

  void write_element_info_f(int *,
			    int *,
			    int *,
			    int *,
			    double *,
			    double *,
			    int *,
			    int *,
			    char *,int);

  void write_global_info_f(char *,
			   int *,
			   int *,
			   int *,
			   int *,
			   int *,
			   int *,
			   int *,
			   char *,
			   char *,
			   int,
			   int,
			   int);

  void write_props_f(double *,
		     double *,
		     int *,
		     int *,
		     int *,
		     int *,
		     int *,
		     int *,
		     int *,
		     int *,
		     char *,
		     char *,
		     int *,
		     char *,
		     int,
		     int,
		     int);

  void write_sparse_info_f(int *,
			   int *,
			   int *,
			   int *,
			   double *,
			   int *,
			   int *,
			   char *,int);

  void write_strscomp_f(double *,
			double *,
			double *,
			int *,
			int *,
			char *,int);

  void write_subiter_f(double *,
		       double *,
		       int *,
		       int *,
		       int *,
		       int *,
		       char *,int);

  void write_ucd_mesh_f(double *,
		       	int *,
		       	int *,
		       	int *,
		       	int *,
		       	int *,
		       	double *,
		       	int *,
		       	int *,
		       	int *,
		       	char *,int);
}

#endif // lithomop3d_externs_h


// version
// $Id: lithomop3d_externs.h,v 1.14 2005/03/12 02:03:18 willic3 Exp $

// End of file
