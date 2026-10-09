#pragma once
#include "Common.h"
#include "Matrix.h"

typedef enum {
	Hartley,		//最优
	Dev,			//稍逊
	Bounding_Box,	//最差
	None			//不归一化
}Normalize_Method;

template<typename _T> struct Point_2D {
	unsigned int m_iCamera_Index;
	unsigned int m_iPoint_Index;
	_T m_Pos[2];
};

template<typename _T>void Get_K4_by_eq_focal(_T eq_f, _T w, _T h, _T K[4]);
template<typename _T>void Get_K3_by_eq_focal(_T eq_f, _T w, _T h, _T K[3]);
template<typename _T>void Get_K9_by_eq_focal(_T eq_f, _T w, _T h, _T K[3 * 3]);
template<typename _T>void K3_Proj(_T K[3], _T P[3], _T uv[2], int bHas_z = 1);
template<typename _T>void K3_Inv(_T K[3], _T K_Inv[3]);	//K3 快速求逆
template<typename _T>void K9_2_K4(_T K9[9], _T K4[4]);
template<typename _T>void K4_2_K9(_T K4[9], _T K9[4]);
template<typename _T>void K3_2_K9(_T K3[9], _T K9[4]);

//一组归一化函数
template<typename _T>void Normalize_by_B_Box_2d(_T P[][2], int n, _T fBox_Size, _T Norm[][2], _T K[3] = NULL, _T K_Inv[] = NULL);	//求一bounding_Box, 将图片搞里头
template<typename _T>void Normalize_Hartley_2d(_T P[][2], int n, _T Norm[][2], _T K[3] = NULL, _T K_Inv[] = NULL);
template<typename _T>void Normalize_Dev_2d(_T P[][2], int n, _T Norm[][2], _T K[4] = NULL, _T K_Inv[4] = NULL);
template<typename _T>void Normalize_2d(_T P[][2], int n, _T Norm[][2],
	Normalize_Method iMethod = Hartley, _T K[4] = NULL, _T K_Inv[4] = NULL);

//H矩阵估计函数
template<typename _T>int Estimate_H_Ref(_T P[][2], _T uv[][2], int n, _T H[3 * 3], int bNormalize = 1, int bUse_SVD = 1);
template<typename _T>int Estimate_H_Zhang(_T Norm_P[][2], _T uv[][2], int n, _T H[3 * 3],
	_T K_Norm[], int bNormalize = 1, Normalize_Method iMethod = Dev);
template<typename _T>void Gen_H_Coeff_z_0(_T P[][2], _T uv[][2], int n, _T A[]);
template<typename _T>_T Test_H_2d(_T P[][2], _T uv[][2], int n, _T H[3 * 3]);

//四元数,旋转矩阵，旋转向量互换函数
template<typename _T>void Quaternion_2_Rotation_Matrix(_T Q[4], _T R[]);	//四元数转旋转矩阵
template<typename _T>void Quaternion_2_Rotation_Vector(_T Q[4], _T V[4]);	//四元数转旋转向量
void Quaternion_Add(float Q_1[], float Q_2[], float Q_3[]);	//四元数加
void Quaternion_Minus(float Q_1[], float Q_2[], float Q_3[]);	//四元数减，其实这两个函数可以用普通向量加减
void Quaternion_Conj(float Q_1[], float Q_2[]);				//简单求个共轭
void Quaternion_Multiply(float Q_1[], float Q_2[], float Q_3[]);	//四元数乘法
void Quaternion_Inv(float Q_1[], float Q_2[]);	//四元数求逆
template<typename _T>void Rotation_Matrix_2_Quaternion(_T R[], _T Q[]);		//旋转矩阵转换四元数
template<typename _T>void Rotation_Matrix_2_Vector(_T R[3 * 3], _T V[4]);		//旋转矩阵转旋转向量
template<typename _T>void Rotation_Matrix_2_Vector_3(_T R[3 * 3], _T V[3]);
template<typename _T>void Rotation_Matrix_2_Vector_4(_T R[3 * 3], _T V[4]);
template<typename _T>void Rotation_Vector_4_2_Matrix(_T V[4], _T R[3 * 3]);		//旋转向量转换旋转矩阵
template<typename _T>void Rotation_Vector_3_2_Matrix(_T V[3], _T R[3 * 3]);		//旋转向量转换旋转矩阵
template<typename _T>void Rotation_Vector_2_Quaternion(_T V[4], _T Q[4]);	//旋转向量转换四元数
template<typename _T>void Rotation_Vector_3_2_4(_T V[], _T V1[]);
template<typename _T>void Rotation_Vector_4_2_3(_T V[], _T V1[]);

//各种投影函数
template<typename _T>int  Get_Pn_by_Pd(_T Pd[2], _T Distort[5], _T Pn[2]);
template<typename _T>void Get_Pd_by_Pn(_T Pn[2], _T D[5], _T Pd);
template<typename _T>void Get_uv_Ref(_T P[4], _T T[4 * 4], _T K[3 * 3], _T D[5], _T uv[2]);
//对于归一化平面的(x,y), 求畸变增量。Pd = Pn + delta, 以下就是求Delta
template<typename _T>void Get_Distort_Coeff(_T x, _T y, _T D[5], _T dc[2 * 5]);
template<typename _T>void Get_Distort_Value(_T x, _T y, _T D[5], _T d[2]);
template<typename _T>void TP(_T T[3 * 4], _T P[3], _T Pt[3]);

//李群李代数
template<typename _T>void Get_dPn_dP(_T Pt[3], _T R[3 * 3], _T dPn_dP[2 * 3]);
template<typename _T>void Get_dPn_dPt(_T Pt[3], _T R[3 * 3], _T dPn_dPt[2 * 3]);	//投影雅可比矩阵
template<typename _T>void Get_dPd_dPn(_T x, _T y, _T k1, _T k2, _T k3, _T p1, _T p2, _T dPd_dPn[4]);
template<typename _T>void Get_PnP_Deriv(_T P[4], _T uv_Ref[2], _T T[4 * 4], _T K[3 * 3], _T D[5],
	_T dE_dK[2*4]=NULL, _T dE_dKsi[2*6] = NULL, _T dE_dD[2*5] = NULL, _T dE_dP[2*3] = NULL, _T E[2]=NULL);
template<typename _T>void Get_dRP_dphi(_T Rp[3], _T dRP_dphi[3 * 3]);
template<typename _T>void Get_dTP_dKsi(_T Pt[3], _T Deriv[4 * 6]);
template<typename _T>void Get_dTP_dKsi(_T T[3 * 4], _T P[2], _T Deriv[4 * 6]);
template<typename _T>void Get_J_E_uv(_T P[4], _T uv[2], _T T[3 * 4], _T K[4], _T D[5], _T J[2][15], _T E[2]);
template<typename _T>void Get_J_E_Norm(_T P[4], _T uv[2], _T T[3 * 4], _T K[4], _T D[5], _T J[2][15], _T E[2]);
template<typename _T>void Gen_Pose_By_R_t(_T R[], _T t[], _T T[],int b4x4 = 1);
template<typename _T>void Get_Inv_T(_T T[3 * 4], _T T_Inv[3 * 4], int b4x4 = 1);
template<typename _T>void T_2_R9_t(_T T[3 * 4], _T R[3 * 3], _T t[3]);
template<typename _T>void Gen_Pose_By_V3_t(_T V3[], _T t[], _T T[]);
template<typename _T>void Get_R_t(_T T[4 * 4], _T R[3 * 3], _T t[3]);
template<typename _T>void Hat(_T V[], _T M[]);
template<typename _T>void Exp_V3(_T V[], _T B[], _T eps=1e-5);
template<typename _T>int bIs_Hat(_T A[], int n);	//判断是否反对称
template<typename _T>void Vee(_T M[], _T V[3]);
template<typename _T>void se3_2_SE3(_T Ksi[6], _T T[]);

//一组BCH 求雅可比函数
template<typename _T>void Get_Jl_3(_T Rotation_Vector[3], _T J[]);
template<typename _T>void Get_Jl_4(_T Rotation_Vector_4[4], _T J[]);
template<typename _T>void Get_Jr_3(_T phi[3], _T J[]);
template<typename _T>void Get_Jr_4(_T phi[4], _T J[]);

//一气呵成构造H 矩阵
template<typename _T>void Get_H_JtE(_T T[][3 * 4], _T K[4], _T D[5], _T P[][2], Point_2D<_T> uv[],
	int iObservation_Count, int iOrder, _T H[], _T JtE[]);
template<typename _T>void Get_H_JtE(_T T[][3 * 4], _T K[4], _T D[5], _T P[][2], Point_2D<_T> uv[],
	int iObservation_Count, int iOrder, _T H[], _T JtE[]);

//求E矩阵
template<typename _T>int Estimate_E_Rt_5Point(_T Point_A[5][2], _T Point_B[5][2], _T R[3 * 3], _T t[3], _T E[3 * 3] = NULL);
template<typename _T>int Estimate_E_T_5Point(_T Point_A[5][2], _T Point_B[5][2], _T T[3 * 4], _T E[3 * 3] = NULL);
template<typename _T>int Estimate_E_8Point(_T Point_A[8][2], _T Point_B[8][2], _T T[3 * 4], _T E[3 * 3]);
template<typename _T>int Estimate_F_7Point(_T Point_A[7][2], _T Point_B[7][2], _T F[3][3 * 3], int* piCount=NULL);
template<typename _T>int Estimate_F_8Point(_T Point_A[8][2], _T Point_B[8][2], _T F[3 * 3]);
template<typename _T>int Estimate_H_4Point(_T Point_A[4][2], _T Point_B[4][2], _T T[3 * 4]=NULL, _T H[3 * 3]=NULL, _T *pfError=NULL);
template<typename _T>int Estimate_E_T_nPoint(_T Point_A[5][2], _T Point_B[5][2], int iCount, int* piInlier, _T T[3 * 4] = NULL, _T E[3 * 3] = NULL);
template<typename _T>int Elect_R_t(_T Point_A[5][2], _T Point_B[5][2], int iPoint_Count, _T R[][3 * 3], _T t[][3], int iRt_Count, int* piBest_E_Index, _T* pfError=NULL);
template<typename _T>int Elect_E(_T Point_A[5][2], _T Point_B[5][2], _T E[][3 * 3], int iCount, _T Best_T[3 * 4], int iSample_Count, int* piBest_E_Index);
template<typename _T>_T Sampson(_T x1[2], _T x2[2], _T E[3 * 3]);	//求Sampson距离
template<typename _T>int Ransac_E(_T Point_A[][2], _T Point_B[][2], int iCount, _T T[3 * 4], _T E[3 * 3], int bUse_5_Point = 1, int bPoint_In_Place = 1);
template<typename _T>int Ransac_H(_T Point_A[][2], _T Point_B[][2], int iCount,	_T T[3 * 4], _T H[3 * 3], int bPoint_In_Place = 0, _T eps = 0.0004);

template<typename _T>int Get_Inlier_Count(_T Point_A[][2], _T Point_B[][2], int iCount, _T T[3 * 4], int bSwap_Inlier_Forward = 1, int iMethod = 2,
	_T* pfError = NULL, _T eps = 0.005208333333333333f * 0.005208333333333333);
template<typename _T>void Sample_XY(_T Point_A[][2], _T Point_B[][2], int iPoint_Count, int iSample_Count);
//int Ransac_Remain_Count(int iInlier, int iSample_Count, int iMin_Sample_Count);
//template<typename _T>void Gen_H_Coeff_row(_T P1[2], _T P2[2], _T row[2 * 9]);
//template<typename _T>void H_2_R_t(_T H[3 * 3], _T R[4][3 * 3], _T t[4][3], _T n[4][3] = NULL);

//三角化
template<typename _T>int Triangulate_Cramer(_T Pn1[3], _T Pn2[3], _T T[4 * 4], _T P[3]);
template<typename _T>int Triangulate_Gauss(_T Pn1[3], _T Pn2[3], _T T1[4 * 4], _T T2[4 * 4], _T P[3]);
template<typename _T>int Triangulate_DLT(_T Pn1[3], _T Pn2[3], _T T1[3 * 4], _T T2[3 * 4], _T P[3]);
template<typename _T> void Test_E(_T E[3 * 3], _T P1[5][2], _T P2[5][2], int iCount);
template<typename _T>_T Test_Triangulate(_T P[3], _T T[4 * 4], _T Point_A[2], _T Point_B[2]);
template<typename _T>_T Test_Triangulate(_T T[4 * 4], _T Point_A[][2], _T Point_B[][2], int iCount);