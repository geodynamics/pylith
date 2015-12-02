/* ---------------------------------------------------------------------- */
/* f1 entry function for isotropic linear elasticity in 3-D.
 *
 * Solution fields = [disp(dim), vel(dim)]
 * Auxiliary fields = [lambda(1), mu(1), initialstress(dim*dim),
 *                     initialstrain(dim*dim)]
 */
void
pylith_fekernels_f1_IsotropicLinearElasticity3D(const PylithInt dim,
						const PylithInt numS,
						const PylithInt numA,
						const PylithInt sOff[],
						const PylithInt sOff_x[],
						const PylithScalar s[],
						const PylithScalar s_t[],
						const PylithScalar s_x[],
						const PylithInt aOff[],
						const PylithInt aOff_x[],
						const PylithScalar a[],
						const PylithScalar a_t[],
						const PylithScalar a_x[],
						const PylithReal t,
						const PylithScalar x[],
						PylithScalar f1[])
{ /* f1_IsotropicLinearElasticity3D */
  const PylithInt _dim = 3;

  const PylithInt _numS = 2;
  const PylithInt i_disp = 0;

  const PylithInt _numA = 4;
  const PylithInt i_lambda = 0;
  const PylithInt i_mu = 1;
  const PylithInt i_istress = 2;
  const PylithInt i_istrain = 3;

  const PylithInt numAVol = 3;
  const PylithInt aOffVol[3] = { aOff[i_lambda], aOff[i_istress],
				 aOff[i_istrain] };
  const PylithInt aOffVol_x[3] = { aOff_x[i_lambda], aOff_x[i_istress],
				   aOff_x[i_istrain] };

  const PylithInt numADev = 3;
  const PylithInt aOffDev[3] = { aOff[i_mu], aOff[i_istress], aOff[i_istrain] };
  const PylithInt aOffDev_x[3] = { aOff_x[i_mu], aOff_x[i_istress],
				   aOff_x[i_istrain] };

  assert(_dim == dim);
  assert(_numS == numS);
  assert(_numA == numA);
  assert(sOff);
  assert(sOff_x);
  assert(aOff);
  assert(aOff_x);

  pylith_fekernels_volumetricStress_IsotropicLinearElasticity3D(dim, 1, numAVol, &sOff[i_disp], &sOff_x[i_disp], s, s_t, s_x, aOffVol, aOffVol_x, a, a_t, a_x, t, x, f1);
  pylith_fekernels_deviatoricStress_IsotropicLinearElasticity3D(dim, 1, numADev, &sOff[i_disp], &sOff_x[i_disp], s, s_t, s_x, aOffDev, aOffDev_x, a, a_t, a_x, t, x, f1);
} /* f1_IsotropicLinearElasticity3D */


/* ---------------------------------------------------------------------- */
/* g3_uu entry function for isotropic linear elasticity in 3-D.
 *
 * Solution fields = [disp(dim), vel(dim)]
 * Auxiliary fields = [lambda(1), mu(1)]
 */
void
pylith_fekernels_g3_uu_IsotropicLinearElasticity3D(const PylithInt dim,
						   const PylithInt numS,
						   const PylithInt numA,
						   const PylithInt sOff[],
						   const PylithInt sOff_x[],
						   const PylithScalar s[],
						   const PylithScalar s_t[],
						   const PylithScalar s_x[],
						   const PylithInt aOff[],
						   const PylithInt aOff_x[],
						   const PylithScalar a[],
						   const PylithScalar a_t[],
						   const PylithScalar a_x[],
						   const PylithReal t,
						   const PylithReal utshift,
						   const PylithScalar x[],
						   PylithScalar g3[])
{ /* g3_uu_IsotropicLinearElasticity3D */
  const PylithInt _dim = 3;

  const PylithInt _numS = 2;

  const PylithInt _numA = 2;
  const PylithInt i_lambda = 0;
  const PylithInt i_mu = 1;

  const PylithScalar lambda = a[aOff[i_lambda]];
  const PylithScalar mu = a[aOff[i_mu]];

  const PylithScalar mu2 = 2.0 * mu;
  const PylithScalar lambda2mu = lambda + mu2;
   
  const PylithReal C1111 = lambda2mu;
  const PylithReal C2222 = lambda2mu;
  const PylithReal C3333 = lambda2mu;
  const PylithReal C1122 = lambda;
  const PylithReal C1133 = lambda;
  const PylithReal C2233 = lambda;
  const PylithReal C1212 = mu2;
  const PylithReal C2323 = mu2;
  const PylithReal C1313 = mu2;

  assert(_dim == dim);
  assert(_numS == numS);
  assert(_numA == numA);
  assert(aOff);

  /*
    g(i,j,di,dj) = C(i,di,j,dj)

     0: g0000 = C1111
     1: g0001 = C1112
     2: g0002 = C1113
     3: g0010 = C1211 symmetry C1112
     4: g0011 = C1212
     5: g0012 = C1213
     6: g0020 = C1311 symmetry C1113
     7: g0021 = C1312 symmetry C1213
     8: g0022 = C1313
  
     9: g0100 = C1121 symmetry C1112
    10: g0101 = C1122
    11: g0102 = C1123
    12: g0110 = C1221 symmetry C1212
    13: g0111 = C1222
    14: g0112 = C1223
    15: g0120 = C1321 symmetry C1213
    16: g0121 = C1322
    17: g0122 = C1323
  
    18: g0200 = C1131 symmetry C1113
    19: g0201 = C1132 symmetry C1123
    20: g0202 = C1133
    21: g0210 = C1231 symmetry C1213
    22: g0211 = C1232 symmetry C1223
    23: g0212 = C1233
    24: g0220 = C1331 symmetry C1313
    25: g0221 = C1332 symmetry C1323
    26: g0222 = C1333
  
    27: g1000 = C2111 symmetry C1112
    28: g1001 = C2112 symmetry C1212
    29: g1002 = C2113 symmetry C1213
    30: g1010 = C2211 symmetry C1122
    31: g1011 = C2212 symmetry C1222
    32: g1012 = C2213 symmetry C1322
    33: g1020 = C2311 symmetry C1123
    34: g1021 = C2312 symmetry C1223
    35: g1022 = C2313 symmetry C1323
  
    36: g1100 = C2121 symmetry C1212
    37: g1101 = C2122 symmetry C1222
    38: g1102 = C2123 symmetry C1223
    39: g1110 = C2221 symmetry C1222
    40: g1111 = C2222
    41: g1112 = C2223
    42: g1120 = C2321 symmetry C1223
    43: g1121 = C2322 symmetry C2223
    44: g1122 = C2323
  
    45: g1200 = C2131 symmetry C1213
    46: g1201 = C2132 symmetry C1223
    47: g1202 = C2133 symmetry C1233
    48: g1210 = C2231 symmetry C1322
    49: g1211 = C2232 symmetry C2223
    50: g1212 = C2233
    51: g1220 = C2331 symmetry C1323
    52: g1221 = C2332 symmetry C2323
    53: g1222 = C2333
  
    54: g2000 = C3111 symmetry C1113
    55: g2001 = C3112 symmetry C1213
    56: g2002 = C3113 symmetry C1313
    57: g2010 = C3211 symmetry C1123
    58: g2011 = C3212 symmetry C1223
    59: g2012 = C3213 symmetry C1323
    60: g2020 = C3311 symmetry C1133
    61: g2021 = C3312 symmetry C1233
    62: g2022 = C3313 symmetry C1333
 
    63: g2100 = C3121 symmetry C1213
    64: g2101 = C3122 symmetry C1322
    65: g2102 = C3123 symmetry C1323
    66: g2110 = C3221 symmetry C1223
    67: g2111 = C3222 symmetry C2223
    68: g2112 = C3223 symmetry C2323
    69: g2120 = C3321 symmetry C1233
    70: g2121 = C3322 symmetry C2233
    71: g2122 = C3323 symmetry C2333
  
    72: g2200 = C3131 symmetry C1313
    73: g2201 = C3132 symmetry C1323
    74: g2202 = C3133 symmetry C1333
    75: g2210 = C3231 symmetry C1323
    76: g2211 = C3232 symmetry C2323
    77: g2212 = C3233 symmetry C2333
    78: g2220 = C3331 symmetry C1333
    79: g2221 = C3332 symmetry C2333
    80: g2222 = C3333
  
  */

  g3[ 0] += C1111; /* g0000 */
  g3[ 4] += C1212; /* g0011 */
  g3[ 8] += C1313; /* g0022 */
  g3[10] += C1122; /* g0101 */
  g3[12] += C1212; /* g0110, C1221 */
  g3[20] += C1133; /* g0202 */
  g3[24] += C1313; /* g0220, C1331 */
  g3[28] += C1212; /* g1001, C2112 */
  g3[30] += C1122; /* g1010, C2211 */
  g3[40] += C2222; /* g1111 */
  g3[44] += C2323; /* g1122 */
  g3[50] += C2233; /* g1212 */
  g3[52] += C2323; /* g1221, C2332 */
  g3[56] += C1313; /* g2002, C3113 */
  g3[60] += C1133; /* g2020, C3311 */
  g3[68] += C2323; /* g2112, C3223 */
  g3[70] += C2233; /* g2121, C3322 */
  g3[80] += C3333; /* g2222 */
} /* g3_uu_IsotropicLinearElasticity3D */


/* ---------------------------------------------------------------------- */
/* Calculate volumetric stress for isotropic linear elasticity.
 *
 * Solution fields = [disp(dim)]
 * Auxiliary fields = [lambda(1), initialstress(dim*dim), initialstrain(dim*dim)]
 *
 * We compute the stress relative to a reference stress/strain state.
 *
 * stress_ij - initialstress_ij = lambda * (strain_kk - initialstrain_kk) * delta_ij + 2*mu * (strain_ij - initialstrain_ij)
 */
void
pylith_fekernels_volumetricStress_IsotropicLinearElasticity3D(const PylithInt dim,
							      const PylithInt numS,
							      const PylithInt numA,
							      const PylithInt sOff[],
							      const PylithInt sOff_x[],
							      const PylithScalar s[],
							      const PylithScalar s_t[],
							      const PylithScalar s_x[],
							      const PylithInt aOff[],
							      const PylithInt aOff_x[],
							      const PylithScalar a[],
							      const PylithScalar a_t[],
							      const PylithScalar a_x[],
							      const PylithReal t,
							      const PylithScalar x[],
							      PylithScalar stress[])
{ /* volumetricStress_IsotropicLinearElasticity3D */
  const PylithInt _dim = 3;

  const PylithInt _numS = 1;
  const PylithInt i_disp = 0;
  const PylithScalar* disp_x = &s_x[sOff[i_disp]];

  const PylithInt _numA = 3;
  const PylithInt i_lambda = 0;
  const PylithInt i_istress = 1;
  const PylithInt i_istrain = 2;
  const PylithScalar lambda = a[aOff[i_lambda]];
  const PylithScalar* initialstress = &a[aOff[i_istress]];
  const PylithScalar* initialstrain = &a[aOff[i_istrain]];

  PylithInt i;
  PylithScalar trace = 0;
  PylithScalar meanistress = 0;

  assert(_dim == dim);
  assert(_numS == numS);
  assert(_numA == numA);
  assert(sOff);
  assert(aOff);
  assert(s_x);
  assert(a);
  assert(stress);

  for (i=0; i < _dim; ++i) {
    trace += disp_x[i] - initialstrain[i*_dim+i];
    meanistress += initialstress[i*_dim+i];
  } /* for */
  meanistress /= (PylithScalar)_dim;
  for (i = 0; i < _dim; ++i) {
    stress[i*_dim+i] += lambda * trace + meanistress;
  } /* for */
} /* volumetricStress_IsotropicLinearElasticity */


/* ---------------------------------------------------------------------- */
/* Calculate deviatoric stress for isotropic linear elasticity.
 *
 * Solution fields = [disp(dim)]
 * Auxiliary fields = [mu(1), initialstress(dim*dim), initialstrain(dim*dim)]
 */
void
pylith_fekernels_deviatoricStress_IsotropicLinearElasticity3D(const PylithInt dim,
							      const PylithInt numS,
							      const PylithInt numA,
							      const PylithInt sOff[],
							      const PylithInt sOff_x[],
							      const PylithScalar s[],
							      const PylithScalar s_t[],
							      const PylithScalar s_x[],
							      const PylithInt aOff[],
							      const PylithInt aOff_x[],
							      const PylithScalar a[],
							      const PylithScalar a_t[],
							      const PylithScalar a_x[],
							      const PylithReal t,
							      const PylithScalar x[],
							      PylithScalar stress[])
{ /* deviatoricStress_IsotropicLinearElasticity3D */
  const PylithInt _dim = 3;

  const PylithInt _numS = 1;
  const PylithInt i_disp = 0;
  const PylithScalar* disp_x = &s_x[sOff[i_disp]];

  const PylithInt _numA = 3;
  const PylithInt i_mu = 0;
  const PylithInt i_istress = 1;
  const PylithInt i_istrain = 2;
  const PylithScalar mu = a[aOff[i_mu]];
  const PylithScalar* initialstress = &a[aOff[i_istress]];
  const PylithScalar* initialstrain = &a[aOff[i_istrain]];

  PylithInt i, j;
  PylithScalar meanistress = 0;

  assert(_dim == dim);
  assert(_numS == numS);
  assert(_numA == numA);
  assert(sOff);
  assert(aOff);
  assert(s_x);
  assert(a);
  assert(stress);

  for (i=0; i < _dim; ++i) {
    meanistress += initialstress[i*_dim+i];
  } /* for */
  meanistress /= (PylithScalar)_dim;
  for (i=0; i < _dim; ++i) {
    for (j=0; j < _dim; ++j) {
      stress[i*_dim+j] += mu * (disp_x[i*_dim+j] + disp_x[j*_dim+i] - initialstrain[i*_dim+j]) + initialstress[i*_dim+j];
    } /* for */
    stress[i*_dim+i] -= meanistress;
  } /* for */
} /* deviatoricStress_IsotropicLinearElasticity3D */


/* ---------------------------------------------------------------------- */
/* f1 entry function for 3-D incompressible isotropic linear elasticity.
 *
 * Solution fields = [disp(dim), pres]
 * Auxiliary fields = [lambda(1), mu(1), initialstress(dim*dim),
 *                     initialstrain(dim*dim)]
 */
void
pylith_fekernels_f1_IncomprUIntegral3D(const PylithInt dim,
				       const PylithInt numS,
				       const PylithInt numA,
				       const PylithInt sOff[],
				       const PylithInt sOff_x[],
				       const PylithScalar s[],
				       const PylithScalar s_t[],
				       const PylithScalar s_x[],
				       const PylithInt aOff[],
				       const PylithInt aOff_x[],
				       const PylithScalar a[],
				       const PylithScalar a_t[],
				       const PylithScalar a_x[],
				       const PylithReal t,
				       const PylithScalar x[],
				       PylithScalar f1[])
{ /* f1_IncomprUIntegral3D */
  const PylithInt _dim = 3;

  const PylithInt _numS = 2;
  const PylithInt i_disp = 0;
  const PylithInt i_pres = 1;
  const PylithScalar pres = s[sOff[i_pres]];

  const PylithInt _numA = 4;

  const PylithInt i_mu = 1;
  const PylithInt i_istress = 2;
  const PylithInt i_istrain = 3;

  const PylithInt numADev = 3;
  const PylithInt aOffDev[3] = { aOff[i_mu], aOff[i_istress], aOff[i_istrain] };
  const PylithInt aOffDev_x[3] = { aOff_x[i_mu], aOff_x[i_istress],
				   aOff_x[i_istrain] };

  assert(_dim == dim);
  assert(_numS == numS);
  assert(_numA == numA);
  assert(sOff);
  assert(sOff_x);
  assert(aOff);
  assert(aOff_x);

  // NOTE:  At present, only the deviatoric initial strains and stresses are
  // being used. Not sure at present how to incorporate the volumetric parts.
  pylith_fekernels_deviatoricStress_IsotropicLinearElasticityIncompr3D(
				dim, 1, numADev, &sOff[i_disp], &sOff_x[i_disp],
				s, s_t, s_x, aOffDev, aOffDev_x, a, a_t, a_x, t,
				x, f1);
  PylithInt i;

  for (i=0; i < _dim; ++i) {
    f1[i*_dim+i] -= pres;
  } /* for */
} /* f1_IncomprUIntegral3D */

/* ---------------------------------------------------------------------- */
/* Calculate deviatoric stress for 3-D incompressible isotropic linear elasticity.
 *
 * Solution fields = [disp(dim)]
 * Auxiliary fields = [mu(1), initialstress(dim*dim), initialstrain(dim*dim)]
 */
void
pylith_fekernels_deviatoricStress_IsotropicLinearElasticityIncompr3D(
					const PylithInt dim,
					const PylithInt numS,
					const PylithInt numA,
					const PylithInt sOff[],
					const PylithInt sOff_x[],
					const PylithScalar s[],
					const PylithScalar s_t[],
					const PylithScalar s_x[],
					const PylithInt aOff[],
					const PylithInt aOff_x[],
					const PylithScalar a[],
					const PylithScalar a_t[],
					const PylithScalar a_x[],
					const PylithReal t,
					const PylithScalar x[],
					PylithScalar stress[])
{ /* deviatoricStress_IsotropicLinearElasticityIncompr3D */
  const PylithInt _dim = 3;

  const PylithInt _numS = 1;
  const PylithInt i_disp = 0;
  const PylithScalar* disp_x = &s_x[sOff[i_disp]];

  const PylithInt _numA = 3;
  const PylithInt i_mu = 0;
  const PylithInt i_istress = 1;
  const PylithInt i_istrain = 2;
  const PylithScalar mu = a[aOff[i_mu]];
  const PylithScalar* initialstress = &a[aOff[i_istress]];
  const PylithScalar* initialstrain = &a[aOff[i_istrain]];

  PylithInt i, j;
  PylithScalar meanistress = 0;
  PylithScalar meanistrain = 0;
  PylithScalar meanstrain = 0;
  PylithScalar devistrain[_dim*_dim];
  PylithScalar devstrain[_dim*_dim];

  assert(_dim == dim);
  assert(_numS == numS);
  assert(_numA == numA);
  assert(sOff);
  assert(aOff);
  assert(s_x);
  assert(a);
  assert(stress);

  // Need to make sure whether all the decomposition into deviatoric parts is
  // necessary.
  for (i=0; i < _dim; ++i) {
    meanistress += initialstress[i*_dim+i];
    meanistrain += initialstrain[i*_dim+i];
    meanstrain += disp_x[i];
  } /* for */
  meanistress /= (PylithScalar)_dim;
  meanistrain /= (PylithScalar)_dim;
  meanstrain /= (PylithScalar)_dim;
  for (i=0; i < _dim; ++i) {
    for (j=0; j < _dim; ++j) {
      devistrain[i*_dim+j] = initialstrain[i*_dim+j];
      devstrain[i*_dim+j] = 0.5 * (disp_x[i*_dim+j] + disp_x[j*dim+i]);
    } /* for */
    devistrain[i*_dim+i] -= meanistrain;
    devstrain[i*_dim+i] -= meanstrain;
  } /* for */
  
  for (i=0; i < _dim; ++i) {
    for (j=0; j < _dim; ++j) {
      stress[i*_dim+j] += mu * (devstrain[i*_dim+j] - devistrain[i*_dim+j]) +
	initialstress[i*_dim+j];
    } /* for */
    stress[i*_dim+i] -= meanistress;
  } /* for */
} /* deviatoricStress_IsotropicLinearElasticityIncompr3D */


/* ---------------------------------------------------------------------- */
/* g3_uu entry function for 3-D incompressible isotropic linear elasticity.
 *
 * Solution fields = [disp(dim), vel(dim)]
 * Auxiliary fields = [lambda(1), mu(1), initialstress(dim*dim),
 *                     initialstrain(dim*dim)]
 */
void
pylith_fekernels_g3_uu_IncompressIsotropicLinearElasticity3D(
					const PylithInt dim,
					const PylithInt numS,
					const PylithInt numA,
					const PylithInt sOff[],
					const PylithInt sOff_x[],
					const PylithScalar s[],
					const PylithScalar s_t[],
					const PylithScalar s_x[],
					const PylithInt aOff[],
					const PylithInt aOff_x[],
					const PylithScalar a[],
					const PylithScalar a_t[],
					const PylithScalar a_x[],
					const PylithReal t,
					const PylithReal utshift,
					const PylithScalar x[],
					PylithScalar g3[])
{ /* g3_uu_IncomprIsotropicLinearElasticity3D */
  const PylithInt _dim = 3;

  const PylithInt _numS = 2;

  const PylithInt _numA = 4;
  const PylithInt i_mu = 1;

  const PylithScalar mu = a[aOff[i_mu]];

  const PylithReal C1111 = 5.0 * mu/3.0;
  const PylithReal C2222 = C1111;
  const PylithReal C3333 = C1111;
  const PylithReal C1212 = mu;
  const PylithReal C1122 = -mu/3.0;

  assert(_dim == dim);
  assert(_numS == numS);
  assert(_numA == numA);
  assert(aOff);


  /* g(f,g,df,dg) = C_{f,df,g,dg}^{\prime} -
                    \frac{1}{6} \sum_h C_{f,df,h,h}^{\prime} \delta_{g,dg}

    00: g0000 = C1111 = 5*mu/3
    01: g0001 = C1112
    02: g0002 = C1113
    03: g0010 = C1211
    04: g0011 = C1212 = mu
    05: g0012 = C1213
    06: g0020 = C1311
    07: g0021 = C1312
    08: g0022 = C1313 = mu
    09: g0100 = C1121
    10: g0101 = C1122 = -mu/3
    11: g0102 = C1123
    12: g0110 = C1221 = mu
    13: g0111 = C1222
    14: g0112 = C1223
    15: g0120 = C1321
    16: g0121 = C1322
    17: g0122 = C1323
    18: g0200 = C1131
    19: g0201 = C1132
    20: g0202 = C1133 = -mu/3
    21: g0210 = C1231
    22: g0211 = C1232
    23: g0212 = C1233
    24: g0220 = C1331 = mu
    25: g0221 = C1332
    26: g0222 = C1333
    27: g1000 = C2111
    28: g1001 = C2112 = mu
    29: g1002 = C2113
    30: g1010 = C2211 = -mu/3
    31: g1011 = C2212
    32: g1012 = C2213
    33: g1020 = C2311
    34: g1021 = C2312
    35: g1022 = C2313
    36: g1100 = C2121 = mu
    37: g1101 = C2122
    38: g1102 = C2123
    39: g1110 = C2221
    40: g1111 = C2222 = 5*mu/3
    41: g1112 = C2223
    42: g1120 = C2321
    43: g1121 = C2322
    44: g1122 = C2323 = mu
    45: g1200 = C2131
    46: g1201 = C2132
    47: g1202 = C2133
    48: g1210 = C2231
    49: g1211 = C2232
    50: g1212 = C2233 = -mu/3
    51: g1220 = C2331
    52: g1221 = C2332 = mu
    53: g1222 = C2333
    54: g2000 = C3111
    55: g2001 = C3112
    56: g2002 = C3113 = mu
    57: g2010 = C3211
    58: g2011 = C3212
    59: g2012 = C3213
    60: g2020 = C3311 = -mu/3
    61: g2021 = C3312
    62: g2022 = C3313
    63: g2100 = C3121
    64: g2101 = C3122
    65: g2102 = C3123
    66: g2110 = C3221
    67: g2111 = C3222
    68: g2112 = C3223 = mu
    69: g2120 = C3321
    70: g2121 = C3322 = -mu/3
    71: g2122 = C3323
    72: g2200 = C3131 = mu
    73: g2201 = C3132
    74: g2202 = C3133
    75: g2210 = C3231
    76: g2211 = C3232 = mu
    77: g2212 = C3233
    78: g2220 = C3331
    79: g2221 = C3332
    80: g2222 = C3333 = 5*mu/3
  */

  g3[ 0] += C1111; /* g0000 */
  g3[ 4] += C1212; /* g0011 */
  g3[ 8] += C1212; /* g0022 */
  g3[10] += C1122; /* g0101, C1221 */
  g3[ 9] += C1212; /* g1001, C2112 */
  g3[10] += C1122; /* g1010, C2211 */
  g3[12] += C1212; /* g1100, C2121 */
  g3[20] += C1122; /* g0202, C1133 */
  g3[24] += C1212; /* g0220, C1331 */
  g3[28] += C1212; /* g1001, C2112 */
  g3[30] += C1122; /* g1010, C2211 */
  g3[36] += C1212; /* g1100, C2121 */
  g3[40] += C2222; /* g1111, C2222 */
  g3[44] += C1212; /* g1122, C2323 */
  g3[50] += C1122; /* g1212, C2233 */
  g3[52] += C1212; /* g1221, C2332 */
  g3[56] += C1212; /* g2002, C3113 */
  g3[60] += C1122; /* g2020, C3311 */
  g3[68] += C1212; /* g2112, C3223 */
  g3[70] += C1122; /* g2121, C3322 */
  g3[72] += C1212; /* g2200, C3131 */
  g3[76] += C1212; /* g2211, C3232 */
  g3[80] += C3333; /* g2222, C3333 */

} /* g3_uu_IncomprIsotropicLinearElasticity3D */
					      

/* ---------------------------------------------------------------------- */
/* g2_uv entry function for 3-D incompressible isotropic linear elasticity.
 *
 * Solution fields = [disp(dim), pres]
 * Auxiliary fields = None
 */
void
pylith_fekernels_g2_uv_IncomprIsotropicLinearElasticity3D(
					const PylithInt dim,
					const PylithInt numS,
					const PylithInt numA,
					const PylithInt sOff[],
					const PylithInt sOff_x[],
					const PylithScalar s[],
					const PylithScalar s_t[],
					const PylithScalar s_x[],
					const PylithInt aOff[],
					const PylithInt aOff_x[],
					const PylithScalar a[],
					const PylithScalar a_t[],
					const PylithScalar a_x[],
					const PylithReal t,
					const PylithReal utshift,
					const PylithScalar x[],
					PylithScalar g2[])
{ /* g2_uv_IncomprIsotropicLinearElasticityPlaneStrain */
  const PylithInt _numS = 2;

  const PylithInt _numA = 0;

  PylithInt i;

  assert(_numS == numS);
  assert(_numA == numA);

  for (i=0; i < dim; ++i) {
    g2[i*dim+i] += +1.0;
  } /* for */
} /* g2_uv_IncomprIsotropicLinearElasticity3D */

