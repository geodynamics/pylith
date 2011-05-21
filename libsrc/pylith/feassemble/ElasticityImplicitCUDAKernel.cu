#include <stdexcept>

#define CHECK_CUDA_ERROR_MSG(e, msg) do {if (e) return e;} while (0)

__global__ void integrateElasticity(float *elemMat, float *geometry, float *analytic)
{
  const int        gridIdx = blockIdx.x + blockIdx.y*gridDim.x; /* Indexes element batch */
  const int        Kidx    = threadIdx.x + threadIdx.y*12; /* This is (i,j) for test and basis functions */
  const int        idx     = Kidx;                        /* Unique thread ID (K block is for a single element) */

  const int        Goffset = gridIdx*288;
  __shared__ float G[288];
  const int        Koffset = Kidx*9;
  float            K[9];
  const int        Eoffset = gridIdx*4608;
  float            E       = 0.0;

  // Load geometry from global memory into G in shared memory
G[idx+0] = geometry[Goffset+idx+0];
G[idx+144] = geometry[Goffset+idx+144];

  /* Copy K^{ij} into local memory (not coalesced) */
  K[0] = analytic[Koffset+0];
  K[1] = analytic[Koffset+1];
  K[2] = analytic[Koffset+2];
  K[3] = analytic[Koffset+3];
  K[4] = analytic[Koffset+4];
  K[5] = analytic[Koffset+5];
  K[6] = analytic[Koffset+6];
  K[7] = analytic[Koffset+7];
  K[8] = analytic[Koffset+8];
  __syncthreads(); /* Make G available */
  /* Do contraction */ 
  /*   NEED TO INTERLEAVE CONTRACTIONS OF CONCURRENT ELEMENTS? See Volkov talk */
  E += G[0] * K[0];
E += G[1] * K[1];
E += G[2] * K[2];
E += G[3] * K[3];
E += G[4] * K[4];
E += G[5] * K[5];
E += G[6] * K[6];
E += G[7] * K[7];
E += G[8] * K[8];

  /* Store contraction result */
  elemMat[Eoffset+idx+0] = E;
  E = 0.0;
  E += G[9] * K[0];
E += G[10] * K[1];
E += G[11] * K[2];
E += G[12] * K[3];
E += G[13] * K[4];
E += G[14] * K[5];
E += G[15] * K[6];
E += G[16] * K[7];
E += G[17] * K[8];

  /* Store contraction result */
  elemMat[Eoffset+idx+144] = E;
  E = 0.0;
  E += G[18] * K[0];
E += G[19] * K[1];
E += G[20] * K[2];
E += G[21] * K[3];
E += G[22] * K[4];
E += G[23] * K[5];
E += G[24] * K[6];
E += G[25] * K[7];
E += G[26] * K[8];

  /* Store contraction result */
  elemMat[Eoffset+idx+288] = E;
  E = 0.0;
  E += G[27] * K[0];
E += G[28] * K[1];
E += G[29] * K[2];
E += G[30] * K[3];
E += G[31] * K[4];
E += G[32] * K[5];
E += G[33] * K[6];
E += G[34] * K[7];
E += G[35] * K[8];

  /* Store contraction result */
  elemMat[Eoffset+idx+432] = E;
  E = 0.0;
  E += G[36] * K[0];
E += G[37] * K[1];
E += G[38] * K[2];
E += G[39] * K[3];
E += G[40] * K[4];
E += G[41] * K[5];
E += G[42] * K[6];
E += G[43] * K[7];
E += G[44] * K[8];

  /* Store contraction result */
  elemMat[Eoffset+idx+576] = E;
  E = 0.0;
  E += G[45] * K[0];
E += G[46] * K[1];
E += G[47] * K[2];
E += G[48] * K[3];
E += G[49] * K[4];
E += G[50] * K[5];
E += G[51] * K[6];
E += G[52] * K[7];
E += G[53] * K[8];

  /* Store contraction result */
  elemMat[Eoffset+idx+720] = E;
  E = 0.0;
  E += G[54] * K[0];
E += G[55] * K[1];
E += G[56] * K[2];
E += G[57] * K[3];
E += G[58] * K[4];
E += G[59] * K[5];
E += G[60] * K[6];
E += G[61] * K[7];
E += G[62] * K[8];

  /* Store contraction result */
  elemMat[Eoffset+idx+864] = E;
  E = 0.0;
  E += G[63] * K[0];
E += G[64] * K[1];
E += G[65] * K[2];
E += G[66] * K[3];
E += G[67] * K[4];
E += G[68] * K[5];
E += G[69] * K[6];
E += G[70] * K[7];
E += G[71] * K[8];

  /* Store contraction result */
  elemMat[Eoffset+idx+1008] = E;
  E = 0.0;
  E += G[72] * K[0];
E += G[73] * K[1];
E += G[74] * K[2];
E += G[75] * K[3];
E += G[76] * K[4];
E += G[77] * K[5];
E += G[78] * K[6];
E += G[79] * K[7];
E += G[80] * K[8];

  /* Store contraction result */
  elemMat[Eoffset+idx+1152] = E;
  E = 0.0;
  E += G[81] * K[0];
E += G[82] * K[1];
E += G[83] * K[2];
E += G[84] * K[3];
E += G[85] * K[4];
E += G[86] * K[5];
E += G[87] * K[6];
E += G[88] * K[7];
E += G[89] * K[8];

  /* Store contraction result */
  elemMat[Eoffset+idx+1296] = E;
  E = 0.0;
  E += G[90] * K[0];
E += G[91] * K[1];
E += G[92] * K[2];
E += G[93] * K[3];
E += G[94] * K[4];
E += G[95] * K[5];
E += G[96] * K[6];
E += G[97] * K[7];
E += G[98] * K[8];

  /* Store contraction result */
  elemMat[Eoffset+idx+1440] = E;
  E = 0.0;
  E += G[99] * K[0];
E += G[100] * K[1];
E += G[101] * K[2];
E += G[102] * K[3];
E += G[103] * K[4];
E += G[104] * K[5];
E += G[105] * K[6];
E += G[106] * K[7];
E += G[107] * K[8];

  /* Store contraction result */
  elemMat[Eoffset+idx+1584] = E;
  E = 0.0;
  E += G[108] * K[0];
E += G[109] * K[1];
E += G[110] * K[2];
E += G[111] * K[3];
E += G[112] * K[4];
E += G[113] * K[5];
E += G[114] * K[6];
E += G[115] * K[7];
E += G[116] * K[8];

  /* Store contraction result */
  elemMat[Eoffset+idx+1728] = E;
  E = 0.0;
  E += G[117] * K[0];
E += G[118] * K[1];
E += G[119] * K[2];
E += G[120] * K[3];
E += G[121] * K[4];
E += G[122] * K[5];
E += G[123] * K[6];
E += G[124] * K[7];
E += G[125] * K[8];

  /* Store contraction result */
  elemMat[Eoffset+idx+1872] = E;
  E = 0.0;
  E += G[126] * K[0];
E += G[127] * K[1];
E += G[128] * K[2];
E += G[129] * K[3];
E += G[130] * K[4];
E += G[131] * K[5];
E += G[132] * K[6];
E += G[133] * K[7];
E += G[134] * K[8];

  /* Store contraction result */
  elemMat[Eoffset+idx+2016] = E;
  E = 0.0;
  E += G[135] * K[0];
E += G[136] * K[1];
E += G[137] * K[2];
E += G[138] * K[3];
E += G[139] * K[4];
E += G[140] * K[5];
E += G[141] * K[6];
E += G[142] * K[7];
E += G[143] * K[8];

  /* Store contraction result */
  elemMat[Eoffset+idx+2160] = E;
  E = 0.0;
  E += G[144] * K[0];
E += G[145] * K[1];
E += G[146] * K[2];
E += G[147] * K[3];
E += G[148] * K[4];
E += G[149] * K[5];
E += G[150] * K[6];
E += G[151] * K[7];
E += G[152] * K[8];

  /* Store contraction result */
  elemMat[Eoffset+idx+2304] = E;
  E = 0.0;
  E += G[153] * K[0];
E += G[154] * K[1];
E += G[155] * K[2];
E += G[156] * K[3];
E += G[157] * K[4];
E += G[158] * K[5];
E += G[159] * K[6];
E += G[160] * K[7];
E += G[161] * K[8];

  /* Store contraction result */
  elemMat[Eoffset+idx+2448] = E;
  E = 0.0;
  E += G[162] * K[0];
E += G[163] * K[1];
E += G[164] * K[2];
E += G[165] * K[3];
E += G[166] * K[4];
E += G[167] * K[5];
E += G[168] * K[6];
E += G[169] * K[7];
E += G[170] * K[8];

  /* Store contraction result */
  elemMat[Eoffset+idx+2592] = E;
  E = 0.0;
  E += G[171] * K[0];
E += G[172] * K[1];
E += G[173] * K[2];
E += G[174] * K[3];
E += G[175] * K[4];
E += G[176] * K[5];
E += G[177] * K[6];
E += G[178] * K[7];
E += G[179] * K[8];

  /* Store contraction result */
  elemMat[Eoffset+idx+2736] = E;
  E = 0.0;
  E += G[180] * K[0];
E += G[181] * K[1];
E += G[182] * K[2];
E += G[183] * K[3];
E += G[184] * K[4];
E += G[185] * K[5];
E += G[186] * K[6];
E += G[187] * K[7];
E += G[188] * K[8];

  /* Store contraction result */
  elemMat[Eoffset+idx+2880] = E;
  E = 0.0;
  E += G[189] * K[0];
E += G[190] * K[1];
E += G[191] * K[2];
E += G[192] * K[3];
E += G[193] * K[4];
E += G[194] * K[5];
E += G[195] * K[6];
E += G[196] * K[7];
E += G[197] * K[8];

  /* Store contraction result */
  elemMat[Eoffset+idx+3024] = E;
  E = 0.0;
  E += G[198] * K[0];
E += G[199] * K[1];
E += G[200] * K[2];
E += G[201] * K[3];
E += G[202] * K[4];
E += G[203] * K[5];
E += G[204] * K[6];
E += G[205] * K[7];
E += G[206] * K[8];

  /* Store contraction result */
  elemMat[Eoffset+idx+3168] = E;
  E = 0.0;
  E += G[207] * K[0];
E += G[208] * K[1];
E += G[209] * K[2];
E += G[210] * K[3];
E += G[211] * K[4];
E += G[212] * K[5];
E += G[213] * K[6];
E += G[214] * K[7];
E += G[215] * K[8];

  /* Store contraction result */
  elemMat[Eoffset+idx+3312] = E;
  E = 0.0;
  E += G[216] * K[0];
E += G[217] * K[1];
E += G[218] * K[2];
E += G[219] * K[3];
E += G[220] * K[4];
E += G[221] * K[5];
E += G[222] * K[6];
E += G[223] * K[7];
E += G[224] * K[8];

  /* Store contraction result */
  elemMat[Eoffset+idx+3456] = E;
  E = 0.0;
  E += G[225] * K[0];
E += G[226] * K[1];
E += G[227] * K[2];
E += G[228] * K[3];
E += G[229] * K[4];
E += G[230] * K[5];
E += G[231] * K[6];
E += G[232] * K[7];
E += G[233] * K[8];

  /* Store contraction result */
  elemMat[Eoffset+idx+3600] = E;
  E = 0.0;
  E += G[234] * K[0];
E += G[235] * K[1];
E += G[236] * K[2];
E += G[237] * K[3];
E += G[238] * K[4];
E += G[239] * K[5];
E += G[240] * K[6];
E += G[241] * K[7];
E += G[242] * K[8];

  /* Store contraction result */
  elemMat[Eoffset+idx+3744] = E;
  E = 0.0;
  E += G[243] * K[0];
E += G[244] * K[1];
E += G[245] * K[2];
E += G[246] * K[3];
E += G[247] * K[4];
E += G[248] * K[5];
E += G[249] * K[6];
E += G[250] * K[7];
E += G[251] * K[8];

  /* Store contraction result */
  elemMat[Eoffset+idx+3888] = E;
  E = 0.0;
  E += G[252] * K[0];
E += G[253] * K[1];
E += G[254] * K[2];
E += G[255] * K[3];
E += G[256] * K[4];
E += G[257] * K[5];
E += G[258] * K[6];
E += G[259] * K[7];
E += G[260] * K[8];

  /* Store contraction result */
  elemMat[Eoffset+idx+4032] = E;
  E = 0.0;
  E += G[261] * K[0];
E += G[262] * K[1];
E += G[263] * K[2];
E += G[264] * K[3];
E += G[265] * K[4];
E += G[266] * K[5];
E += G[267] * K[6];
E += G[268] * K[7];
E += G[269] * K[8];

  /* Store contraction result */
  elemMat[Eoffset+idx+4176] = E;
  E = 0.0;
  E += G[270] * K[0];
E += G[271] * K[1];
E += G[272] * K[2];
E += G[273] * K[3];
E += G[274] * K[4];
E += G[275] * K[5];
E += G[276] * K[6];
E += G[277] * K[7];
E += G[278] * K[8];

  /* Store contraction result */
  elemMat[Eoffset+idx+4320] = E;
  E = 0.0;
  E += G[279] * K[0];
E += G[280] * K[1];
E += G[281] * K[2];
E += G[282] * K[3];
E += G[283] * K[4];
E += G[284] * K[5];
E += G[285] * K[6];
E += G[286] * K[7];
E += G[287] * K[8];

  /* Store contraction result */
  elemMat[Eoffset+idx+4464] = E;
  
  
}

cudaError_t setupKernel(const int dim, const int numBasisFuncs, const int numCells, float *K, float **geometry, float **elemMat, float **analytic_gpu, float **geometry_gpu, float **elemMat_gpu)
{
  const int   N             = numCells;
  const int   numComponents = dim;
  int         Ksize         = (numBasisFuncs*numComponents * dim)*(numBasisFuncs*numComponents * dim);
  size_t      Kbytes        = Ksize * sizeof(float);
  int         Gsize         = N*dim*dim;
  size_t      Gbytes        = Gsize * sizeof(float);
  int         Esize         = (numBasisFuncs*numComponents)*(numBasisFuncs*numComponents);
  size_t      Ebytes        = Esize * sizeof(float);
  cudaError_t cerr;

  cerr = cudaMalloc(analytic_gpu, Kbytes);CHECK_CUDA_ERROR_MSG(cerr, "CUDA Allocation failure");
  cerr = cudaMemcpy(*analytic_gpu, K, Kbytes, cudaMemcpyHostToDevice);CHECK_CUDA_ERROR_MSG(cerr, "CUDA Allocation failure");
  cerr = cudaMallocHost(geometry, Gbytes, cudaHostAllocDefault);CHECK_CUDA_ERROR_MSG(cerr, "CUDA Allocation failure");
  cerr = cudaMalloc(geometry_gpu, Gbytes);CHECK_CUDA_ERROR_MSG(cerr, "CUDA Allocation failure");
  cerr = cudaMallocHost(elemMat, Ebytes, cudaHostAllocDefault);CHECK_CUDA_ERROR_MSG(cerr, "CUDA Allocation failure");
  cerr = cudaMalloc(elemMat_gpu, Ebytes);CHECK_CUDA_ERROR_MSG(cerr, "CUDA Allocation failure");
  for(int i = 0; i < Esize; ++i) {(*elemMat)[i] = 0.0;}
  return cudaSuccess;
}

// Calculate a conforming thread grid for N kernels
void calculateGrid(const int N, const int blockSize, unsigned int& x, unsigned int& y, unsigned int& z)
{
  z = 1;
  if (N % blockSize) {
    // 'Invalid block size '+str(blockSize)+' for '+str(N)+' elements'
    throw std::runtime_error("Invalid block size");
  }
  const int Nblocks = N/blockSize;
  for(x = (int) (sqrt(Nblocks) + 0.5); x > 0; --x) {
    y = Nblocks/x;
    if (x*y == Nblocks) break;
  }
  if (x*y != Nblocks) {
    // 'Could not find partition for '+str(N)+' with block size '+str(blockSize)
    throw std::runtime_error("Could not find partition");
  }
  return;
}

cudaError_t launchKernel(const int spaceDim, const int numBasis, const int elementBatchSize, const int numConcurrentElements,
                         const int N, float *analytic_gpu, float *geometry_gpu, float *elemMat_gpu)
{
  dim3 grid, block;
  block.x = numBasis*spaceDim;
  block.y = numBasis*spaceDim;
  block.z = numConcurrentElements;
  calculateGrid(N, elementBatchSize, grid.x, grid.y, grid.z);
  // self.logPrint('Running %d elements with Thread Block size %d' % (N, reduce(int.__mul__, blockDim)), debugLevel=1, debugSection=self.section)
  // self.logPrint('  using grid '+str(self.calculateGrid(N, elementBatchSize)), debugLevel=1, debugSection=self.section)
  integrateElasticity<<<grid, block>>>(elemMat_gpu, geometry_gpu, analytic_gpu);
  return cudaSuccess;
}

cudaError_t cleanupKernel(float *geometry, float *elemMat, float *analytic_gpu, float *geometry_gpu, float *elemMat_gpu)
{
  cudaError_t cerr;

  cerr = cudaFree(analytic_gpu);CHECK_CUDA_ERROR_MSG(cerr, "CUDA Deallocation failure");
  cerr = cudaFree(geometry_gpu);CHECK_CUDA_ERROR_MSG(cerr, "CUDA Deallocation failure");
  cerr = cudaFree(elemMat_gpu);CHECK_CUDA_ERROR_MSG(cerr, "CUDA Deallocation failure");
  cerr = cudaFreeHost(geometry);CHECK_CUDA_ERROR_MSG(cerr, "CUDA Deallocation failure");
  cerr = cudaFreeHost(elemMat);CHECK_CUDA_ERROR_MSG(cerr, "CUDA Deallocation failure");
  return cudaSuccess;
}
