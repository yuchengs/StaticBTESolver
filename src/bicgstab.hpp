#include <cassert>
#include <cuda_runtime.h>
#include <cusparse_v2.h>
#include <cublas_v2.h>
#include <iostream>

#define THREAD_ID threadIdx.x+blockIdx.x*blockDim.x
#define THREAD_COUNT gridDim.x*blockDim.x

__global__ void pre_spmv(int n, const double *precondition, const double *vec, double *result) {
    for (unsigned int i = THREAD_ID; i < n; i+=THREAD_COUNT)
        result[i] = precondition[i] * vec[i];
}
__global__ void gpu_VVequalrhat(int num, double *result, double * inval){
    for(unsigned int i = THREAD_ID; i < num; i += THREAD_COUNT) {
        result[i] = inval[i] * 6.8 + 0.6;
    }
}

void print(double* dev_ptr, int sz) {
    double* host_ptr;
    host_ptr = (double*)malloc(sizeof(double) * sz);
    cudaMemcpy(host_ptr, dev_ptr, sz * sizeof(double), cudaMemcpyDeviceToHost);
    for (int i = 0; i < sz; i++) {
        std::cout << host_ptr[i] << " ";
    }
    std::cout << std::endl;
    if (host_ptr) free(host_ptr);
}

class BICGSTAB {
private:
    int dim, max_iters;
    double tol;

    // CUDA sdk context
    cusparseHandle_t cusparse_handle = nullptr;
    cublasHandle_t cublas_handle = nullptr;

    // CUDA sdk status
    cusparseStatus_t cusparse_latest_status = CUSPARSE_STATUS_SUCCESS;
    cublasStatus_t cublas_latest_status = CUBLAS_STATUS_SUCCESS;
    cudaError_t cuda_latest_status = cudaSuccess;

    // temp data
    cusparseMatDescr_t descrA = nullptr;
    int* dev_A_csrRowPtr = nullptr;
    int* dev_A_csrColInd = nullptr;
    double* dev_A_csrVal = nullptr;
    double* dev_b = nullptr;
    int A_nnz;
    double* dev_jacobi = nullptr;

public:
    BICGSTAB(int n, int maxit, double tol=1e-6) : dim(n), max_iters(maxit), tol(tol), A_nnz(0) {
        // cusparse
        this->cusparse_latest_status = cusparseCreate(&this->cusparse_handle);
        assert(CUSPARSE_STATUS_SUCCESS == this->cusparse_latest_status);
        this->cusparse_latest_status = cusparseCreateMatDescr(&this->descrA);
        assert(CUSPARSE_STATUS_SUCCESS == this->cusparse_latest_status);
        cusparseSetMatIndexBase(descrA, CUSPARSE_INDEX_BASE_ZERO);
        cusparseSetMatType(descrA, CUSPARSE_MATRIX_TYPE_GENERAL);
        // cublas
        this->cublas_latest_status = cublasCreate(&this->cublas_handle);
        assert(CUBLAS_STATUS_SUCCESS == this->cublas_latest_status);
    }
    void init(int* csrRowPtr, int* csrColInd, double* csrVal, int nnz, double* b) {
        this->A_nnz = nnz;

        if (dev_A_csrRowPtr) cudaFree(dev_A_csrRowPtr);
        if (dev_A_csrColInd) cudaFree(dev_A_csrColInd);
        if (dev_A_csrVal) cudaFree(dev_A_csrVal);
        if (dev_b) cudaFree(dev_b);

        // allocate device memory
        this->cuda_latest_status = cudaMalloc((void**)&dev_A_csrRowPtr, sizeof(int) * (this->dim + 1));
        assert(cudaSuccess == this->cuda_latest_status);
        this->cuda_latest_status = cudaMalloc((void**)&dev_A_csrColInd, sizeof(int) * this->A_nnz);
        assert(cudaSuccess == this->cuda_latest_status);
        this->cuda_latest_status = cudaMalloc((void**)&dev_A_csrVal, sizeof(double) * this->A_nnz);
        assert(cudaSuccess == this->cuda_latest_status);
        this->cuda_latest_status = cudaMalloc((void**)&dev_b, sizeof(double) * dim);
        assert(cudaSuccess == this->cuda_latest_status);
        // copy from host to device
        this->cuda_latest_status = cudaMemcpy(dev_A_csrRowPtr, csrRowPtr, sizeof(int) * (dim + 1), cudaMemcpyHostToDevice);
        assert(cudaSuccess == this->cuda_latest_status);
        this->cuda_latest_status = cudaMemcpy(dev_A_csrColInd, csrColInd, sizeof(int) * this->A_nnz, cudaMemcpyHostToDevice);
        assert(cudaSuccess == this->cuda_latest_status);
        this->cuda_latest_status = cudaMemcpy(dev_A_csrVal, csrVal, sizeof(double) * this->A_nnz, cudaMemcpyHostToDevice);
        assert(cudaSuccess == this->cuda_latest_status);
        this->cuda_latest_status = cudaMemcpy(dev_b, b, sizeof(double) * dim, cudaMemcpyHostToDevice);
        assert(cudaSuccess == this->cuda_latest_status);
        double* jacobi = nullptr;
        jacobi = (double*)malloc(sizeof(double) * this->dim);
        memset(jacobi, 0, sizeof(double) & this->dim);
        for (int row = 0; row < this->dim; row++){
            for(int j = csrRowPtr[row]; j < csrRowPtr[row + 1]; j++) {
                if (csrColInd[j] == row){
                    jacobi[row] = 1.0 / csrVal[j];
                }
            }
        }
        if (this->dev_jacobi) cudaFree(this->dev_jacobi);
        cudaMalloc((void**)&this->dev_jacobi, sizeof(double) * this->dim);
        cudaMemcpy(this->dev_jacobi, jacobi, sizeof(double) * this->dim, cudaMemcpyHostToDevice);
        free(jacobi);
    }

    void solve(double* x) {
        const double one = 1.0;
        const double zero = 0.0;
        // ########################################################################################################
        // #02: bicgstab begins
        // ########################################################################################################
        // ## declare and initialize dev_x and dev_r
        // TODO: should have copied initial x from host
        double* dev_x = nullptr;
        this->cuda_latest_status = cudaMalloc((void**)&dev_x, sizeof(double) * this->dim);
        assert(cudaSuccess == this->cuda_latest_status);
        double* dev_r = nullptr;
        this->cuda_latest_status = cudaMalloc((void**)&dev_r, sizeof(double) * this->dim);
        assert(cudaSuccess == this->cuda_latest_status);
        this->cuda_latest_status = cudaMemset(dev_x, 0, sizeof(double) * this->dim);
        assert(cudaSuccess == this->cuda_latest_status);
        this->cuda_latest_status = cudaMemset(dev_r, 0, sizeof(double) * this->dim);
        assert(cudaSuccess == this->cuda_latest_status);
        double* dev_t = nullptr;
        this->cuda_latest_status = cudaMalloc((void**)&dev_t, sizeof(double) * this->dim);
        assert(cudaSuccess == this->cuda_latest_status);
        double* dev_ph = nullptr;
        this->cuda_latest_status = cudaMalloc((void**)&dev_ph, sizeof(double) * this->dim);
        assert(cudaSuccess == this->cuda_latest_status);
        double* dev_z = nullptr;
        this->cuda_latest_status = cudaMalloc((void**)&dev_z, sizeof(double) * this->dim);
        assert(cudaSuccess == this->cuda_latest_status);

        cublasDaxpy(this->cublas_handle, this->dim, &one, this->dev_b, 1, dev_r, 1);
        // ### 2: p = r and \tilde{r} = r
        double* dev_p = nullptr;
        this->cuda_latest_status = cudaMalloc((void**)&dev_p, sizeof(double) * this->dim);
        assert(cudaSuccess == this->cuda_latest_status);
        cudaMemset(dev_p, 0, sizeof(double) * this->dim);
        double* dev_rw = nullptr;
        this->cuda_latest_status = cudaMalloc((void**)&dev_rw, sizeof(double) * this->dim);
        assert(cudaSuccess == this->cuda_latest_status);
        cublasDcopy(this->cublas_handle, this->dim, dev_r, 1, dev_p, 1);
        cublasDcopy(this->cublas_handle, this->dim, dev_r, 1, dev_rw, 1);
        double nrmr0, nrmr;
        cublasDnrm2(this->cublas_handle, this->dim, dev_r, 1, &nrmr0);
        nrmr = nrmr0;
        // TODO check abs tol
        double* dev_q = nullptr;
        this->cuda_latest_status = cudaMalloc((void**)&dev_q, sizeof(double) * this->dim);
        assert(cudaSuccess == this->cuda_latest_status);
        cudaMemset(dev_q, 0, sizeof(double) * this->dim);

        double rhop = 1, rho = 1, omega = 1.0, alpha = 1.0, neg_alpha, neg_omega;
        double temp, temp2, beta;
        for (int i = 1; i <= this->max_iters; i++) {
            if (nrmr < 1e-50) break;
            // ### 4: \rho = \tilde{r}^T r
            rhop = rho;
            cublasDdot(this->cublas_handle, this->dim, dev_rw, 1, dev_r, 1, &rho);
            if (rho <1e-200 && rho > -1e-200) {
                rho = 1;
                alpha = 1;
                omega = 1;
                cudaMemset(dev_x, 0, sizeof(double) * this->dim);
                cudaMemset(dev_q, 0, sizeof(double) * this->dim);
                cudaMemset(dev_p, 0, sizeof(double) * this->dim);
                cublasDcopy(this->cublas_handle, this->dim, this->dev_b, 1, dev_r, 1);
                gpu_VVequalrhat<<<24, 1024>>>(this->dim, dev_rw, dev_b);
                continue;
            }
            // ### 12: \beta = (\rho_i / \rho_{i-1} (\alpha / \omega)
            beta = (rho / rhop) * (alpha / omega);
            // ### 13: p = r + \beta (p - \omega q)
            neg_omega = -omega;
            cublasDaxpy(this->cublas_handle, this->dim, &neg_omega, dev_q, 1, dev_p, 1);
            cublasDscal(this->cublas_handle, this->dim, &beta, dev_p, 1);
            cublasDaxpy(this->cublas_handle, this->dim, &one, dev_r, 1, dev_p, 1);
            // ### 15: M \hat{p} = p
            pre_spmv<<<24, 1024>>>(this->dim, this->dev_jacobi, dev_p, dev_ph);
            // ### 16: q = A \hat{p}
            size_t bufferSizeInBytes2;
            void* buffer2;
            cusparseCsrmvEx_bufferSize(this->cusparse_handle,
                                       CUSPARSE_ALG_MERGE_PATH,
                                       CUSPARSE_OPERATION_NON_TRANSPOSE,
                                       this->dim, this->dim,
                                       this->A_nnz,
                                       &one, CUDA_R_64F,
                                       this->descrA,
                                       this->dev_A_csrVal, CUDA_R_64F, this->dev_A_csrRowPtr, this->dev_A_csrColInd,
                                       dev_ph, CUDA_R_64F,
                                       &zero, CUDA_R_64F,
                                       dev_q, CUDA_R_64F,
                                       CUDA_R_64F,
                                       &bufferSizeInBytes2);
            this->cuda_latest_status = cudaMalloc((void**)&buffer2, bufferSizeInBytes2);
            assert(cudaSuccess == this->cuda_latest_status);
            cusparseCsrmvEx(this->cusparse_handle,
                            CUSPARSE_ALG_MERGE_PATH,
                            CUSPARSE_OPERATION_NON_TRANSPOSE,
                            this->dim, this->dim,
                            this->A_nnz,
                            &one, CUDA_R_64F,
                            this->descrA,
                            this->dev_A_csrVal, CUDA_R_64F, this->dev_A_csrRowPtr, this->dev_A_csrColInd,
                            dev_ph, CUDA_R_64F,
                            &zero, CUDA_R_64F,
                            dev_q, CUDA_R_64F,
                            CUDA_R_64F,
                            buffer2);
            cudaFree(buffer2);
            // ### 17: \alpha = \rho_i / (\tilde{r}^T q)
            cublasDdot(this->cublas_handle, this->dim, dev_rw, 1, dev_q, 1, &temp);
            if (temp < 1e-200 && temp > -1e-200) {
                rho = 1;
                alpha = 1;
                omega = 1;
                cudaMemset(dev_x, 0, sizeof(double) * this->dim);
                cudaMemset(dev_q, 0, sizeof(double) * this->dim);
                cudaMemset(dev_p, 0, sizeof(double) * this->dim);
                cublasDcopy(this->cublas_handle, this->dim, this->dev_b, 1, dev_r, 1);
                cublasDcopy(this->cublas_handle, this->dim, this->dev_b, 1, dev_rw, 1);
                gpu_VVequalrhat<<<24, 1024>>>(this->dim, dev_rw, dev_b);
                continue;
            }

            alpha = rho / temp;

            // ### 18: s = r - \alpha q
            neg_alpha = -alpha;
            cublasDaxpy(this->cublas_handle, this->dim, &neg_alpha, dev_q, 1, dev_r, 1);

            pre_spmv<<<24, 1024>>>(this->dim, this->dev_jacobi, dev_r, dev_z);

            // ### 24: t = A \hat{s}
            size_t bufferSizeInBytes3;
            void* buffer3;
            cusparseCsrmvEx_bufferSize(this->cusparse_handle,
                                       CUSPARSE_ALG_MERGE_PATH,
                                       CUSPARSE_OPERATION_NON_TRANSPOSE,
                                       this->dim, this->dim,
                                       this->A_nnz,
                                       &one, CUDA_R_64F,
                                       this->descrA,
                                       this->dev_A_csrVal, CUDA_R_64F, this->dev_A_csrRowPtr, this->dev_A_csrColInd,
                                       dev_z, CUDA_R_64F,
                                       &zero, CUDA_R_64F,
                                       dev_t, CUDA_R_64F,
                                       CUDA_R_64F,
                                       &bufferSizeInBytes3);
            this->cuda_latest_status = cudaMalloc((void**)&buffer3, bufferSizeInBytes3);
            assert(cudaSuccess == this->cuda_latest_status);
            cusparseCsrmvEx(this->cusparse_handle,
                            CUSPARSE_ALG_MERGE_PATH,
                            CUSPARSE_OPERATION_NON_TRANSPOSE,
                            this->dim, this->dim,
                            this->A_nnz,
                            &one, CUDA_R_64F,
                            this->descrA,
                            this->dev_A_csrVal, CUDA_R_64F, this->dev_A_csrRowPtr, this->dev_A_csrColInd,
                            dev_z, CUDA_R_64F,
                            &zero, CUDA_R_64F,
                            dev_t, CUDA_R_64F,
                            CUDA_R_64F,
                            buffer3);
            cudaFree(buffer3);

            // ### 25: \omega = (t^T s) / (t^T t)
            cublasDdot(this->cublas_handle, this->dim, dev_t, 1, dev_r, 1, &temp);
            cublasDdot(this->cublas_handle, this->dim, dev_t, 1, dev_t, 1, &temp2);
            omega = temp / temp2;

            // ### 19: x = x + \alpha \hat{p}
            cublasDaxpy(this->cublas_handle, this->dim, &alpha, dev_ph, 1, dev_x, 1);
            // ### 26: x = x + \omega \hat{s}
            cublasDaxpy(this->cublas_handle, this->dim, &omega, dev_z, 1, dev_x, 1);
            // ### 27: r = s -\omega t
            neg_omega = -omega;
            cublasDaxpy(this->cublas_handle, this->dim, &neg_omega, dev_t, 1, dev_r, 1);
            // ### 20: check for convergence
            cublasDnrm2(this->cublas_handle, this->dim, dev_r, 1, &nrmr);
            if (nrmr / nrmr0 < this->tol) {
                break;
            }
        }
        //std::cout << nrmr << " / " << nrmr0 << " = " << nrmr / nrmr0 << std::endl;
        cudaMemcpy(x, dev_x, this->dim * sizeof(double), cudaMemcpyDeviceToHost);
        // clean up
        if (dev_t) cudaFree(dev_t);
        if (dev_z) cudaFree(dev_z);
        if (dev_ph) cudaFree(dev_ph);
        if (dev_q) cudaFree(dev_q);
        if (dev_p) cudaFree(dev_p);
        if (dev_rw) cudaFree(dev_rw);
        if (dev_x) cudaFree(dev_x);
        if (dev_r) cudaFree(dev_r);
    }
    ~BICGSTAB() {
        if (this->dev_jacobi) cudaFree(this->dev_jacobi);
        if (this->dev_A_csrRowPtr) cudaFree(this->dev_A_csrRowPtr);
        if (this->dev_A_csrColInd) cudaFree(this->dev_A_csrColInd);
        if (this->dev_A_csrVal) cudaFree(this->dev_A_csrVal);
        if (this->dev_b) cudaFree(this->dev_b);

        if (this->descrA) cusparseDestroyMatDescr(this->descrA);
        if (this->cusparse_handle) cusparseDestroy(this->cusparse_handle);
        if (this->cublas_handle) cublasDestroy(this->cublas_handle);
        cudaDeviceReset();
    }
//    void solve(double* x) {
//        const double one = 1.0;
//        const double neg_one = -1.0;
//        const double zero = 0.0;
//        // #01: ilu precondition
//        // ## declare and allocate LU memory
//
//        // ## solve for L and U
//        double* dev_M_csrVal = nullptr;
//        this->cuda_latest_status = cudaMalloc((void**)&dev_M_csrVal, sizeof(double) * this->A_nnz);
//        assert(cudaSuccess == this->cuda_latest_status);
//        this->cuda_latest_status = cudaMemcpy(dev_M_csrVal, this->dev_A_csrVal, sizeof(double) * this->A_nnz, cudaMemcpyDeviceToDevice);
//        assert(cudaSuccess == this->cuda_latest_status);
//
//        // set up L U M misc
//        const cusparseSolvePolicy_t policyM = CUSPARSE_SOLVE_POLICY_NO_LEVEL;
//        const cusparseSolvePolicy_t policyL = CUSPARSE_SOLVE_POLICY_NO_LEVEL;
//        const cusparseSolvePolicy_t policyU = CUSPARSE_SOLVE_POLICY_USE_LEVEL;
//        const cusparseOperation_t transL  = CUSPARSE_OPERATION_NON_TRANSPOSE;
//        const cusparseOperation_t transU  = CUSPARSE_OPERATION_NON_TRANSPOSE;
//
//        // set up L U M description
//        cusparseMatDescr_t descrM = nullptr;
//        cusparseMatDescr_t descrL = nullptr;
//        cusparseMatDescr_t descrU = nullptr;
//        cusparseCreateMatDescr(&descrM);
//        cusparseSetMatIndexBase(descrM, CUSPARSE_INDEX_BASE_ZERO);
//        cusparseSetMatType(descrM, CUSPARSE_MATRIX_TYPE_GENERAL);
//        cusparseCreateMatDescr(&descrL);
//        cusparseSetMatIndexBase(descrL, CUSPARSE_INDEX_BASE_ZERO);
//        cusparseSetMatType(descrL, CUSPARSE_MATRIX_TYPE_GENERAL);
//        cusparseSetMatFillMode(descrL, CUSPARSE_FILL_MODE_LOWER);
//        cusparseSetMatDiagType(descrL, CUSPARSE_DIAG_TYPE_UNIT);
//        cusparseCreateMatDescr(&descrU);
//        cusparseSetMatIndexBase(descrU, CUSPARSE_INDEX_BASE_ZERO);
//        cusparseSetMatType(descrU, CUSPARSE_MATRIX_TYPE_GENERAL);
//        cusparseSetMatFillMode(descrU, CUSPARSE_FILL_MODE_UPPER);
//        cusparseSetMatDiagType(descrU, CUSPARSE_DIAG_TYPE_NON_UNIT);
//
//        // set up L U M info
//        csrilu02Info_t infoM  = nullptr;
//        csrsv2Info_t  infoL  = nullptr;
//        csrsv2Info_t  infoU  = nullptr;
//        cusparseCreateCsrilu02Info(&infoM);
//        cusparseCreateCsrsv2Info(&infoL);
//        cusparseCreateCsrsv2Info(&infoU);
//
//        // query workspace size
//        int pBufferSizeM;
//        int pBufferSizeL;
//        int pBufferSizeU;
//        int pBufferSize;
//        void *pBuffer = nullptr;
//        cusparseDcsrilu02_bufferSize(this->cusparse_handle,
//                                     this->dim, this->A_nnz,
//                                     this->descrA,
//                                     dev_M_csrVal, this->dev_A_csrRowPtr, this->dev_A_csrColInd,
//                                     infoM, &pBufferSizeM);
//        cusparseDcsrsv2_bufferSize(this->cusparse_handle, transL,
//                                   this->dim, this->A_nnz,
//                                   descrL,
//                                   dev_M_csrVal, this->dev_A_csrRowPtr, this->dev_A_csrColInd,
//                                   infoL, &pBufferSizeL);
//        cusparseDcsrsv2_bufferSize(this->cusparse_handle, transU,
//                                   this->dim, this->A_nnz,
//                                   descrU,
//                                   dev_M_csrVal, this->dev_A_csrRowPtr, this->dev_A_csrColInd,
//                                   infoU, &pBufferSizeU);
//        pBufferSize = mmax(pBufferSizeM, mmax(pBufferSizeL, pBufferSizeU));
//
//        // allocate workspace buffer
//        // pBuffer returned by cudaMalloc is automatically aligned to 128 bytes.
//        this->cuda_latest_status = cudaMalloc((void**)&pBuffer, pBufferSize);
//        assert(cudaSuccess == this->cuda_latest_status);
//
//        // analyze for ilu and triangular solve
//        cusparseDcsrilu02_analysis(this->cusparse_handle,
//                                   this->dim, this->A_nnz, descrM,
//                                   dev_M_csrVal, this->dev_A_csrRowPtr, this->dev_A_csrColInd,
//                                   infoM, policyM, pBuffer);
//        //int structural_zero;
//        //this->cusparse_latest_status = cusparseXcsrilu02_zeroPivot(this->cusparse_handle, infoM, &structural_zero);
//        //assert(CUSPARSE_STATUS_SUCCESS == this->cusparse_latest_status);
//        cusparseDcsrsv2_analysis(this->cusparse_handle,
//                                 transL, this->dim, this->A_nnz, descrL,
//                                 dev_M_csrVal, this->dev_A_csrRowPtr, this->dev_A_csrColInd,
//                                 infoL, policyL, pBuffer);
//        cusparseDcsrsv2_analysis(this->cusparse_handle,
//                                 transU, this->dim, this->A_nnz, descrU,
//                                 dev_M_csrVal, this->dev_A_csrRowPtr, this->dev_A_csrColInd,
//                                 infoU, policyU, pBuffer);
//
//        // ## ilu
//        cusparseDcsrilu02(this->cusparse_handle,
//                          this->dim, this->A_nnz, descrA,
//                          dev_M_csrVal, this->dev_A_csrRowPtr, this->dev_A_csrColInd,
//                          infoM, policyM, pBuffer);
////        this->cusparse_latest_status = cusparseXcsrilu02_zeroPivot(this->cusparse_handle, infoM, &structural_zero);
//        //assert(CUSPARSE_STATUS_SUCCESS == this->cusparse_latest_status);
//
//        // ########################################################################################################
//        // #02: bicgstab begins
//        // ########################################################################################################
//        // ## declare and initialize dev_x and dev_r
//        // TODO: should have copied initial x from host
//        double* dev_x = nullptr;
//        this->cuda_latest_status = cudaMalloc((void**)&dev_x, sizeof(double) * this->dim);
//        assert(cudaSuccess == this->cuda_latest_status);
//        double* dev_r = nullptr;
//        this->cuda_latest_status = cudaMalloc((void**)&dev_r, sizeof(double) * this->dim);
//        assert(cudaSuccess == this->cuda_latest_status);
//        this->cuda_latest_status = cudaMemset(dev_x, 0, sizeof(double) * this->dim);
//        assert(cudaSuccess == this->cuda_latest_status);
//        this->cuda_latest_status = cudaMemset(dev_r, 0, sizeof(double) * this->dim);
//        assert(cudaSuccess == this->cuda_latest_status);
//        double* dev_t = nullptr;
//        this->cuda_latest_status = cudaMalloc((void**)&dev_t, sizeof(double) * this->dim);
//        assert(cudaSuccess == this->cuda_latest_status);
//
//        // ### 1: r = b - A x
////        size_t bufferSizeInBytes1;
////        void* buffer1;
////        cusparseCsrmvEx_bufferSize(this->cusparse_handle,
////                                   CUSPARSE_ALG_MERGE_PATH,
////                                   CUSPARSE_OPERATION_NON_TRANSPOSE,
////                                   this->dim, this->dim,
////                                   this->A_nnz,
////                                   &one, CUDA_R_64F,
////                                   this->descrA,
////                                   this->dev_A_csrVal, CUDA_R_64F, this->dev_A_csrRowPtr, this->dev_A_csrColInd,
////                                   dev_x, CUDA_R_64F,
////                                   &zero, CUDA_R_64F,
////                                   dev_r, CUDA_R_64F,
////                                   CUDA_R_64F,
////                                   &bufferSizeInBytes1);
////        this->cuda_latest_status = cudaMalloc((void**)&buffer1, bufferSizeInBytes1);
////        assert(cudaSuccess == this->cuda_latest_status);
////        cusparseCsrmvEx(this->cusparse_handle,
////                        CUSPARSE_ALG_MERGE_PATH,
////                        CUSPARSE_OPERATION_NON_TRANSPOSE,
////                        this->dim, this->dim,
////                        this->A_nnz,
////                        &one, CUDA_R_64F,
////                        this->descrA,
////                        this->dev_A_csrVal, CUDA_R_64F, this->dev_A_csrRowPtr, this->dev_A_csrColInd,
////                        dev_x, CUDA_R_64F,
////                        &zero, CUDA_R_64F,
////                        dev_r, CUDA_R_64F,
////                        CUDA_R_64F,
////                        buffer1);
////        cudaFree(buffer1);
////        cublasDscal(this->cublas_handle, this->dim, &neg_one, dev_r, 1);
//        cublasDaxpy(this->cublas_handle, this->dim, &one, this->dev_b, 1, dev_r, 1);
//        // ### 2: p = r and \tilde{r} = r
//        double* dev_p = nullptr;
//        this->cuda_latest_status = cudaMalloc((void**)&dev_p, sizeof(double) * this->dim);
//        assert(cudaSuccess == this->cuda_latest_status);
//        double* dev_rw = nullptr;
//        this->cuda_latest_status = cudaMalloc((void**)&dev_rw, sizeof(double) * this->dim);
//        assert(cudaSuccess == this->cuda_latest_status);
//        cublasDcopy(this->cublas_handle, this->dim, dev_r, 1, dev_p, 1);
//        cublasDcopy(this->cublas_handle, this->dim, dev_r, 1, dev_rw, 1);
//        double nrmr0, nrmr;
//        cublasDnrm2(this->cublas_handle, this->dim, dev_r, 1, &nrmr0);
//        nrmr = nrmr0;
//        // TODO check abs tol
//        double* dev_q = nullptr;
//        this->cuda_latest_status = cudaMalloc((void**)&dev_q, sizeof(double) * this->dim);
//        assert(cudaSuccess == this->cuda_latest_status);
//        double* dev_ph = nullptr;
//        this->cuda_latest_status = cudaMalloc((void**)&dev_ph, sizeof(double) * this->dim);
//        assert(cudaSuccess == this->cuda_latest_status);
//        double* dev_s = nullptr;
//        this->cuda_latest_status = cudaMalloc((void**)&dev_s, sizeof(double) * this->dim);
//        assert(cudaSuccess == this->cuda_latest_status);
//
//        double rhop = 1, rho = 1, omega = 1.0, alpha = 1.0, beta = 0, neg_alpha, neg_omega;
//        double temp, temp2;
//        for (int i = 0; i < this->max_iters; i++) {
//            if (nrmr < 1e-50) break;
//            // ### 4: \rho = \tilde{r}^T r
//            rhop = rho;
//            cublasDdot(this->cublas_handle, this->dim, dev_rw, 1, dev_r, 1, &rho);
//            if (i > 0) {
//                // ### 12: \beta = (\rho_i / \rho_{i-1} (\alpha / \omega)
//                beta = (rho / rhop) * (alpha / omega);
//                // ### 13: p = r + \beta (p - \omega v)
//                neg_omega = -omega;
//                cublasDaxpy(this->cublas_handle, this->dim, &neg_omega, dev_q, 1, dev_p, 1);
//                cublasDscal(this->cublas_handle, this->dim, &beta, dev_p, 1);
//                cublasDaxpy(this->cublas_handle, this->dim, &one, dev_r, 1, dev_p, 1);
//            }
//            // ### 15: M \hat{p} = p
//            cusparseDcsrsv2_solve(this->cusparse_handle,
//                                  transL, this->dim, this->A_nnz,
//                                  &one, descrL,
//                                  dev_M_csrVal, this->dev_A_csrRowPtr, this->dev_A_csrColInd,
//                                  infoL,
//                                  dev_p,
//                                  dev_t,
//                                  policyL,
//                                  pBuffer);
//            cusparseDcsrsv2_solve(this->cusparse_handle,
//                                  transU, this->dim, this->A_nnz,
//                                  &one, descrU,
//                                  dev_M_csrVal, this->dev_A_csrRowPtr, this->dev_A_csrColInd,
//                                  infoU,
//                                  dev_t,
//                                  dev_ph,
//                                  policyU,
//                                  pBuffer);
//            // ### 16: q = A \hat{p}
//            size_t bufferSizeInBytes2;
//            void* buffer2;
//            cusparseCsrmvEx_bufferSize(this->cusparse_handle,
//                                       CUSPARSE_ALG_MERGE_PATH,
//                                       CUSPARSE_OPERATION_NON_TRANSPOSE,
//                                       this->dim, this->dim,
//                                       this->A_nnz,
//                                       &one, CUDA_R_64F,
//                                       this->descrA,
//                                       this->dev_A_csrVal, CUDA_R_64F, this->dev_A_csrRowPtr, this->dev_A_csrColInd,
//                                       dev_ph, CUDA_R_64F,
//                                       &zero, CUDA_R_64F,
//                                       dev_q, CUDA_R_64F,
//                                       CUDA_R_64F,
//                                       &bufferSizeInBytes2);
//            this->cuda_latest_status = cudaMalloc((void**)&buffer2, bufferSizeInBytes2);
//            assert(cudaSuccess == this->cuda_latest_status);
//            cusparseCsrmvEx(this->cusparse_handle,
//                            CUSPARSE_ALG_MERGE_PATH,
//                            CUSPARSE_OPERATION_NON_TRANSPOSE,
//                            this->dim, this->dim,
//                            this->A_nnz,
//                            &one, CUDA_R_64F,
//                            this->descrA,
//                            this->dev_A_csrVal, CUDA_R_64F, this->dev_A_csrRowPtr, this->dev_A_csrColInd,
//                            dev_ph, CUDA_R_64F,
//                            &zero, CUDA_R_64F,
//                            dev_q, CUDA_R_64F,
//                            CUDA_R_64F,
//                            buffer2);
//            cudaFree(buffer2);
//            // ### 17: \alpha = \rho_i / (\tilde{r}^T q)
//            cublasDdot(this->cublas_handle, this->dim, dev_rw, 1, dev_q, 1, &temp);
//            alpha = rho / temp;
//            // ### 18: s = r - \alpha q
//            neg_alpha = -alpha;
//            cublasDaxpy(this->cublas_handle, this->dim, &neg_alpha, dev_q, 1, dev_r, 1);
//            // ### 19: x = x + \alpha \hat{p}
//            cublasDaxpy(this->cublas_handle, this->dim, &alpha, dev_ph, 1, dev_x, 1);
//            // ### 20: check for convergence
//            cublasDnrm2(this->cublas_handle, this->dim, dev_r, 1, &nrmr);
//            if (nrmr / nrmr0 < this->tol) {
//                break;
//            }
//            // ### 23: M \hat{s} = r
//            cusparseDcsrsv2_solve(this->cusparse_handle,
//                                  transL, this->dim, this->A_nnz,
//                                  &one, descrL,
//                                  dev_M_csrVal, this->dev_A_csrRowPtr, this->dev_A_csrColInd,
//                                  infoL,
//                                  dev_r,
//                                  dev_t,
//                                  policyL,
//                                  pBuffer);
//            cusparseDcsrsv2_solve(this->cusparse_handle,
//                                  transU, this->dim, this->A_nnz,
//                                  &one, descrU,
//                                  dev_M_csrVal, this->dev_A_csrRowPtr, this->dev_A_csrColInd,
//                                  infoU,
//                                  dev_t,
//                                  dev_s,
//                                  policyU,
//                                  pBuffer);
//            // ### 24: t = A \hat{s}
//            size_t bufferSizeInBytes3;
//            void* buffer3;
//            cusparseCsrmvEx_bufferSize(this->cusparse_handle,
//                                       CUSPARSE_ALG_MERGE_PATH,
//                                       CUSPARSE_OPERATION_NON_TRANSPOSE,
//                                       this->dim, this->dim,
//                                       this->A_nnz,
//                                       &one, CUDA_R_64F,
//                                       this->descrA,
//                                       this->dev_A_csrVal, CUDA_R_64F, this->dev_A_csrRowPtr, this->dev_A_csrColInd,
//                                       dev_s, CUDA_R_64F,
//                                       &zero, CUDA_R_64F,
//                                       dev_t, CUDA_R_64F,
//                                       CUDA_R_64F,
//                                       &bufferSizeInBytes3);
//            this->cuda_latest_status = cudaMalloc((void**)&buffer3, bufferSizeInBytes3);
//            assert(cudaSuccess == this->cuda_latest_status);
//            cusparseCsrmvEx(this->cusparse_handle,
//                            CUSPARSE_ALG_MERGE_PATH,
//                            CUSPARSE_OPERATION_NON_TRANSPOSE,
//                            this->dim, this->dim,
//                            this->A_nnz,
//                            &one, CUDA_R_64F,
//                            this->descrA,
//                            this->dev_A_csrVal, CUDA_R_64F, this->dev_A_csrRowPtr, this->dev_A_csrColInd,
//                            dev_s, CUDA_R_64F,
//                            &zero, CUDA_R_64F,
//                            dev_t, CUDA_R_64F,
//                            CUDA_R_64F,
//                            buffer3);
//            cudaFree(buffer3);
//            // ### 25: \omega = (t^T s) / (t^T t)
//            cublasDdot(this->cublas_handle, this->dim, dev_t, 1, dev_r, 1, &temp);
//            cublasDdot(this->cublas_handle, this->dim, dev_t, 1, dev_t, 1, &temp2);
//            omega = temp / temp2;
//            // ### 26: x = x + \omega \hat{s}
//            cublasDaxpy(this->cublas_handle, this->dim, &omega, dev_s, 1, dev_x, 1);
//            // ### 27: r = s -\omega t
//            neg_omega = -omega;
//            cublasDaxpy(this->cublas_handle, this->dim, &neg_omega, dev_t, 1, dev_r, 1);
//            // check for convergence
//            cublasDnrm2(this->cublas_handle, this->dim, dev_r, 1, &nrmr);
//            if (nrmr / nrmr0 < tol) {
//                break;
//            }
//        }
//        std::cout << nrmr << " / " << nrmr0 << " = " << nrmr / nrmr0 << std::endl;
//        cudaMemcpy(x, dev_x, this->dim * sizeof(double), cudaMemcpyDeviceToHost);
//        // clean up
//        if (dev_t) cudaFree(dev_t);
//        if (dev_s) cudaFree(dev_s);
//        if (dev_ph) cudaFree(dev_ph);
//        if (dev_q) cudaFree(dev_q);
//        if (dev_p) cudaFree(dev_p);
//        if (dev_rw) cudaFree(dev_rw);
//        if (dev_x) cudaFree(dev_x);
//        if (dev_r) cudaFree(dev_r);
//        cusparseDestroyMatDescr(descrM);
//        cusparseDestroyMatDescr(descrL);
//        cusparseDestroyMatDescr(descrU);
//        cusparseDestroyCsrilu02Info(infoM);
//        cusparseDestroyCsrsv2Info(infoL);
//        cusparseDestroyCsrsv2Info(infoU);
//        if (dev_M_csrVal) cudaFree(dev_M_csrVal);
//        if (pBuffer) cudaFree(pBuffer);
//    }
};
