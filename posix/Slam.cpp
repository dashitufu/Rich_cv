#include "Image.h"
#include "Slam.h"
template void K9_2_K4(float K9[9], float K4[4]);
template void K9_2_K4(double K9[9], double K4[4]);
template<typename _T>void K9_2_K4(_T K9[9],_T K4[4])
{
	K4[0] = K9[0];
	K4[1] = K9[4];
	K4[2] = K9[2];
	K4[3] = K9[5];
}

template void K4_2_K9(float K4[9], float K9[4]);
template void K4_2_K9(double K4[9], double K9[4]);
template<typename _T>void K4_2_K9(_T K4[9], _T K9[4])
{
	memset(K9, 0, 9 * sizeof(_T));
	K9[0] = K4[0];
	K9[4] = K4[1];
	K9[2] = K4[2];
	K9[5] = K4[3];
	K9[8] = 1;
}

template void K3_2_K9(double K3[9], double K9[4]);
template<typename _T>void K3_2_K9(_T K3[9], _T K9[4])
{
	memset(K9, 0, 9 * sizeof(_T));
	K9[0] = K9[4] = K3[0];
	K9[2] = K3[1];
	K9[5] = K3[2];
	K9[8] = 1;
}

template void Get_K9_by_eq_focal(float eq_f, float w, float h, float K[3 * 3]);
template void Get_K9_by_eq_focal(double eq_f, double w, double h, double K[3 * 3]);
template<typename _T>void Get_K9_by_eq_focal(_T eq_f, _T w, _T h, _T K[3 * 3])
{//通过等效焦距求相机内参矩阵K
	//eq_f： equivalent focal length 等效焦距，单位mm，对应长宽36mm * 24mm 成像平面
	//w, h 相机图片的长宽
	//例：Get_K<_T>(35, 4000, 2256, K);	35mm等效小孔成像焦距，4000x2256希昂宿平面
	memset(K, 0, 9 * sizeof(_T));
	K[0] = K[4] = (_T)w * eq_f / 36;
	K[2] = (_T)w / 2;
	K[5] = (_T)h / 2;

	//此处就是变为(u,v,1)齐次坐标的分别所在
	//_T z = w;
	//K[8] = (_T)z * 35 / 36;	//非齐次
	K[8] = 1;				//(u,v)必齐次
	return;
}

//简化版
template void Get_K4_by_eq_focal(double eq_f, double w, double h, double K[4]);
template<typename _T>void Get_K4_by_eq_focal(_T eq_f, _T w, _T h, _T K[4])
{//通过等效焦距求相机内参矩阵K
	//eq_f： equivalent focal length 等效焦距，单位mm，对应长宽36mm * 24mm 成像平面
	//w, h 相机图片的长宽
	//例：Get_K<_T>(35, 4000, 2256, K);	35mm等效小孔成像焦距，4000x2256希昂宿平面
	K[0] = K[1] = (_T)w * eq_f / 36;
	K[2] = (_T)w / 2;
	K[3] = (_T)h / 2;
	return;
}

template void Get_K3_by_eq_focal(double eq_f, double w, double h, double K[3]);
template<typename _T>void Get_K3_by_eq_focal(_T eq_f, _T w, _T h, _T K[3])
{//通过等效焦距求相机内参矩阵K
	//eq_f： equivalent focal length 等效焦距，单位mm，对应长宽36mm * 24mm 成像平面
	//w, h 相机图片的长宽
	//例：Get_K<_T>(35, 4000, 2256, K);	35mm等效小孔成像焦距，4000x2256希昂宿平面
	K[0] = (_T)w * eq_f / 36;
	K[1] = (_T)w / 2;
	K[2] = (_T)h / 2;
	return;
}

template void K3_Proj(float K[3], float P[3], float uv[2], int bHas_z);
template void K3_Proj(double K[3], double P[3], double uv[2], int bHas_z);
template<typename _T>void K3_Proj(_T K[3], _T P[3], _T uv[],int bHas_z)
{//求空间一点P经过相机内参K投影到成像平面的值
//注意：K3 做不到尺度不变性，只是用来简化计算，一乘scale 就打乱
	if (bHas_z)
	{//注意，本来K 是个正经的三维投影到二维的矩阵
		uv[0] = (K[0] * P[0]) / P[2] + K[1];
		uv[1] = (K[0] * P[1]) / P[2] + K[2];
		//= (0 * px + 0 * py + 1*pz)/(0 * px + 0 * py + 1*pz) =1;
		uv[2] = 1;
	}else
	{//干二维的缩放平移只是副业，此时假想P为齐次坐标(x,y,1)
		uv[0] = K[0] * P[0] + K[1];
		uv[1] = K[0] * P[1] + K[2];
	}
	return;
}

template void K3_Inv(float K[3], float K_Inv[3]);
template void K3_Inv(double K[3], double K_Inv[3]);
template<typename _T>void K3_Inv(_T K[3], _T K_Inv[3])
{//求K3投影矩阵的逆
//将 K的最简形式拆开，u = (f * px)/pz + cx
//	px = pz * (u - cx) / f = pz * (u * 1 / f - cx / f)
//	py = pz * (v - cy) / f = pz * (v * 1 / f - cy / f)
//	= >
//	px = pz *	1/f		0		-cx/f	*	u
//	py			0		1/f		-cy/f		v
//	1			0		0		1			1
	K_Inv[0] = 1.f / K[0];
	K_Inv[1] = -K[1] / K[0];
	K_Inv[2] = -K[2] / K[0];
	return;
}

template<typename _T>void K4_Inv(_T K[4], _T K_Inv[4])
{//求K3投影矩阵的逆
//将 K的最简形式拆开，	u = (fx * px)/pz + cx
//						v = (fy * py)/pz + cy
//	px = pz * (u - cx) / fx = pz * (u * 1 / f - cx / f)
//	py = pz * (v - cy) / fy = pz * (v * 1 / f - cy / f)
//	= >
//	px = pz *	1/fx	0		-cx/fx *	u
//	py			0		1/fx	-cy/fy		v
//	1			0		0		1			1
	K_Inv[0] = 1.f / K[0];
	K_Inv[1] = 1.f / K[1];
	K_Inv[2] = -K[2] / K[0];
	K_Inv[3] = -K[3] / K[1];
}

template<typename _T>void Homo_Normalize(_T V0[], int n, _T V1[])
{//homogeneous normalization, 其实就是归一化为齐次坐标，屁营养也美哟iiu
	_T fCoeff = V0[n - 1];
	for (int i = 0; i < n - 1; i++)
		V1[i] = V0[i] / fCoeff;
	return;
}

template void Normalize_by_B_Box_2d(float P[][2], int n, float fBox_Size, float Norm[][2], float K[3], float K_Inv[]);
template void Normalize_by_B_Box_2d(double P[][2], int n, double fBox_Size, double Norm[][2], double K[3], double K_Inv[]);
template<typename _T>void Normalize_by_B_Box_2d(_T P[][2], int n,
	_T fBox_Size, _T Norm[][2], _T K[3], _T K_Inv[])
{//给所有点求个Bouding Box，将中心对准中心，相当于投影到一个平面上
//实践证明，这个方法虽然能降低条件数，但是远比Hartlay 方法差，本身
//bounding box 就是把最差的点算进来，这些烂点扭曲了数据的分布形态
	_T B_Box[2][2];
	Get_Bounding_Box(P, n, B_Box);

	_T Box_Size[2] = { B_Box[1][0] - B_Box[0][0], B_Box[1][1] - B_Box[0][1] };
	_T fMax = Max(Box_Size[0], Box_Size[1]);
	_T Center[2] = { B_Box[0][0] + Box_Size[0] / 2,
					B_Box[0][1] + Box_Size[1] / 2 };

	_T fScale = fBox_Size / fMax;
	// org * s + delta = new => delta = new - org*s 
	_T Offset[2] = { fBox_Size / 2 - Center[0] * fScale, fBox_Size / 2 - Center[1] * fScale };
	_T K1[3] = { fScale,Offset[0],Offset[1] };

	for (int i = 0; i < n; i++)
		K3_Proj<_T>(K1, P[i], Norm[i], 0);

	if (K)
		memcpy(K, K1, 3 * sizeof(_T));
	if (K_Inv)
		K3_Inv(K1, K_Inv);
	return;
}

template void Normalize_Hartley_2d(float P[][2], int n, float Norm[][2], float K[3], float K_Inv[]);
template void Normalize_Hartley_2d(double P[][2], int n, double Norm[][2], double K[3], double K_Inv[]);
template<typename _T>void Normalize_Hartley_2d(_T P[][2], int n,
	_T Norm[][2], _T K[3], _T K_Inv[])
{//用Hartley 归一化方法
	const _T sqrt_2 = (_T)sqrt(2);

	//先找期望
	_T E[2];
	Get_E_2d(P, n, E);

	_T fDev = 0;
	for (int i = 0; i < n; i++)
		fDev += (_T)sqrt(fGet_Distance(P[i], E, 2));
	fDev /= n;

	_T K1[3] = { sqrt_2 / fDev,-K1[0] * E[0], -K1[0] * E[1] };

	for (int i = 0; i < n; i++)
		K3_Proj(K1, P[i], Norm[i], 0);

	if (K)
		memcpy(K, K1, 3 * sizeof(_T));

	if (K_Inv)
		K3_Inv(K1, K_Inv);
	return;
}
template void Normalize_Dev_2d(float P[][2], int n, float Norm[][2], float K[4], float K_Inv[4]);
template void Normalize_Dev_2d(double P[][2], int n, double Norm[][2], double K[4], double K_Inv[4]);
template<typename _T>void Normalize_Dev_2d(_T P[][2], int n, _T Norm[][2], _T K[4], _T K_Inv[4])
{//用偏离绝对值归一化，此乃安全函数，源可以等于目的
	//注意投影概念， K：是uv 投影到 Norm 平面上的投影矩阵
	//K_Inv 是由Norm 平面恢复到uv 上的投影矩阵
	//K， K_Inv 都是压缩存储

	_T E[2];
	Get_E_2d(P, n, E);

	_T Dev[2] = { 0 };
	Get_Dev_2d<_T>(P, n, E, Dev);

	_T s[2] = { 1.f / Dev[0],  1.f / Dev[1] };
	for (int i = 0; i < n; i++)
	{
		Norm[i][0] = s[0] * (P[i][0] - E[0]);
		Norm[i][1] = s[1] * (P[i][1] - E[1]);
	}

	if(K || K_Inv)
	{
		_T K1[4] = { s[0],s[1], -E[0] * s[0], -E[1] * s[1] };
		if(K)
			memcpy(K, K1, 4 * sizeof(_T));
		if (K_Inv)
			K4_Inv(K1, K_Inv);
	}
	return;
}

template<typename _T>static void Normalize_Zhang(_T Point_2D[][2], int iPoint_Count,
	_T Norm_Point[][2], _T Scale[2], _T Offset[2])
{//对一组点归一化，这个和Colmap又不一样，不求Max,Min
	//毫无营养
	_T Mean[2] = { 0 };
	int i;
	//对该点集中求中心
	for (i = 0; i < iPoint_Count; i++)
	{
		//printf("%f %f\n", Point_2D[i][0], Point_2D[i][1]);
		Mean[0] += Point_2D[i][0], Mean[1] += Point_2D[i][1];
	}

	Mean[0] /= iPoint_Count;
	Mean[1] /= iPoint_Count;

	//注意，这不是标准差，叫1：平均绝对偏差；2，L1 Loss；3，Mean Absolute Erro
	_T Dev[2] = {};	//再求个偏离度
	for (i = 0; i < iPoint_Count; i++)
		Dev[0] += abs(Point_2D[i][0] - Mean[0]), Dev[1] += abs(Point_2D[i][1] - Mean[1]);
	Dev[0] /= iPoint_Count;
	Dev[1] /= iPoint_Count;

	//将绝对偏差的导数作为scale，对每个样本求一个偏离度，再乘以这个scale
	_T s[2] = { 1.f / Dev[0],  1.f / Dev[1] };
	for (i = 0; i < iPoint_Count; i++)
	{
		Norm_Point[i][0] = s[0] * (Point_2D[i][0] - Mean[0]);
		Norm_Point[i][1] = s[1] * (Point_2D[i][1] - Mean[1]);
	}
	Scale[0] = s[0], Scale[1] = s[1];
	Offset[0] = -Mean[0] * s[0], Offset[1] = -Mean[1] * s[1];
	return;
}

template void Gen_H_Coeff_z_0(float P[][2], float uv[][2], int n, float A[]);
template void Gen_H_Coeff_z_0(double P[][2], double uv[][2], int n, double A[]);
template<typename _T>void Gen_H_Coeff_z_0(_T P[][2], _T uv[][2], int n, _T A[])
{//很简单，就是给 z=0 的平面点集构造系数矩阵A
	_T* A1 = A, * P_Cur = P[0], * uv_Cur = uv[0];
	for (int i = 0; i < n; i++)
	{
		_T px = P_Cur[0], py = P_Cur[1], u = uv_Cur[0], v = uv_Cur[1];
		//x	y	1	0	0	0	-ux	-uy	-u
		A1[0] = px, A1[1] = py, A1[2] = 1;
		A1[3] = A1[4] = A1[5] = 0;
		A1[6] = -px * u, A1[7] = -py * u, A1[8] = -u;
		A1 += 9;

		//0	0	0	x	y	1	-vx	-vy	-v
		A1[0] = A1[1] = A1[2] = 0;
		A1[3] = px, A1[4] = py, A1[5] = 1;
		A1[6] = -px * v, A1[7] = -py * v, A1[8] = -v;
		A1 += 9;
		P_Cur += 2, uv_Cur += 2;
	}
}
template int Estimate_H_Zhang(float Norm_P[][2], float uv[][2], int n, float H[3 * 3],	float K_Norm[], int bNormalize, Normalize_Method iMethod);
template int Estimate_H_Zhang(double Norm_P[][2], double uv[][2], int n, double H[3 * 3],	double K_Norm[], int bNormalize, Normalize_Method iMethod);
template<typename _T>int Estimate_H_Zhang(_T Norm_P[][2], _T uv[][2], int n, _T H[3 * 3], 
	_T K_Norm[], int bNormalize, Normalize_Method iMethod)
{//尽可能快
	//**********先分配**************************************/
	int iSize = ALIGN_SIZE_8(n * 2 * 9 * sizeof(_T)) +
		ALIGN_SIZE_8(n * 2 * sizeof(_T)) * 2, bRet = 0;
	_T* A = NULL, (*pNorm_uv)[2];

	unsigned char* p = (unsigned char*)pMalloc(iSize);
	if (!p)
		goto END;
	A = (_T*)p, p += ALIGN_SIZE_8(n * 2 * 9 * sizeof(_T));
	pNorm_uv = (_T(*)[2])p;
	//**********先分配**************************************/

	_T K_uv_Inv[4];
	if (bNormalize)
	{
		Normalize_2d<_T>(uv, n, pNorm_uv, iMethod, NULL, K_uv_Inv);
		Gen_H_Coeff_z_0(Norm_P, pNorm_uv, n, A);
	}else
		Gen_H_Coeff_z_0(Norm_P, uv, n, A);
	//Disp((_T*)pNorm_uv, n, 2, "Norm");

	//******************接着解齐次矛盾方程组 Ax =0*****************/
	_T AtA[9 * 9];
	Transpose_Multiply(A, n * 2, 9, AtA, 0);
	//Disp(AtA, 9, 9, "AtA");
	//Disp(A, n * 2, 9);
	int iResult;
	iResult = bInverse_Power<_T>(AtA, 9, (_T*)NULL, H, (_T)1e-8);
	if (!iResult)
		Solve_Homo_Linear_SVD(A, n * 2, 9, H, &iResult);

	if (!iResult)
		goto END;
	//******************接着解齐次矛盾方程组 Ax =0*****************/

	if (bNormalize)
	{//此时需要恢复Scale
		if (iMethod == Normalize_Method::Dev)
		{
			//Normalize by Dev
			_T K_uv_Inv_1[3 * 3] = { K_uv_Inv[0],0,K_uv_Inv[2],
								0,K_uv_Inv[1],K_uv_Inv[3],
								0,0,1 };
			_T K_P_1[3 * 3] = { K_Norm[0],0,K_Norm[2],
								0,K_Norm[1],K_Norm[3],
								0,0,1 };
			Matrix_Multiply_3x3<_T>(K_uv_Inv_1, H, H);
			Matrix_Multiply_3x3<_T>(H, K_P_1, H);
		}else if (iMethod == Normalize_Method::Hartley ||
			iMethod == Normalize_Method::Bounding_Box)
		{
			//Hartley 方法
			_T K_uv_Inv_1[3 * 3] = { K_uv_Inv[0],0,K_uv_Inv[1],
								0,K_uv_Inv[0],K_uv_Inv[2],
								0,0,1 };
			_T K_P_1[3 * 3] = { K_Norm[0],0,K_Norm[1],
								0,K_Norm[0],K_Norm[2],
								0,0,1 };
			Matrix_Multiply_3x3<_T>(K_uv_Inv_1, H, H);
			Matrix_Multiply_3x3<_T>(H, K_P_1, H);
		}else
			printf("Not implemented\n");
	}
	//Disp(H, 3, 3, "H");
	//Normalize(H, 9, H);
	bRet = 1;
END:
	Free(A);
	return bRet;
}

template void Normalize_2d(float P[][2], int n, float Norm[][2], Normalize_Method iMethod, float K[4], float K_Inv[4]);
template void Normalize_2d(double P[][2], int n, double Norm[][2], Normalize_Method iMethod, double K[4], double K_Inv[4]);
template<typename _T>void Normalize_2d(_T P[][2], int n, _T Norm[][2],
	Normalize_Method iMethod, _T K[4], _T K_Inv[4])
{
	switch (iMethod)
	{
	case Dev:
		Normalize_Dev_2d<_T>(P, n, Norm, K, K_Inv);
		break;
	case Hartley:
		Normalize_Hartley_2d(P, n, Norm, K, K_Inv);
		break;
	case Bounding_Box:
		Normalize_by_B_Box_2d<_T>(P, n, 1, Norm, K, K_Inv);
		break;
	case None:
		memcpy(Norm, P, n * 2 * sizeof(_T));
		break;
	default:
		printf("Not implemented in Normalize_2d\n");
	}
	return;
}

template int Estimate_H_Ref(float P[][2], float uv[][2], int n, float H[3 * 3], int bNormalize, int bUse_SVD);
template int Estimate_H_Ref(double P[][2], double uv[][2], int n, double H[3 * 3], int bNormalize, int bUse_SVD);
template<typename _T>int Estimate_H_Ref(_T P[][2], _T uv[][2], int n, _T H[3 * 3], int bNormalize, int bUse_SVD)
{//估计一个H矩阵，满足(u,v,1)' = 1/z * HP
//这是有限制的H矩阵估计，要求参考点Point_Ref 同一定死在z=0 平面上
	const Normalize_Method Norm_Method = Normalize_Method::Dev;

	//**********先分配**************************************/
	int iSize = ALIGN_SIZE_8(n * 2 * 9 * sizeof(_T)) +
		ALIGN_SIZE_8(n * 2 * sizeof(_T)) * 2, bRet = 0;
	_T* A = NULL, (*pNorm_P)[2], (*pNorm_uv)[2];
	unsigned char* p = (unsigned char*)pMalloc(iSize);
	if (!p)
		goto END;

	A = (_T*)p, p += ALIGN_SIZE_8(n * 2 * 9 * sizeof(_T));
	pNorm_P = (_T(*)[2])p;	p += ALIGN_SIZE_8(n * 2 * sizeof(_T));
	pNorm_uv = (_T(*)[2])p;
	//**********先分配**************************************/

	_T K_P[4], K_uv_Inv[4];
	if (bNormalize)
	{
		if(Norm_Method == Normalize_Method::Dev)
		{
			Normalize_Dev_2d<_T>(P, n, pNorm_P, K_P);
			Normalize_Dev_2d<_T>(uv, n, pNorm_uv, NULL, K_uv_Inv);
		}else if(Norm_Method == Normalize_Method::Hartley)
		{
			Normalize_Hartley_2d<_T>(P, n, pNorm_P, K_P);
			Normalize_Hartley_2d<_T>(uv, n, pNorm_uv, NULL, K_uv_Inv);
		}else if(Norm_Method == Normalize_Method::Bounding_Box)
		{
			Normalize_by_B_Box_2d<_T>(P, n,1, pNorm_P, K_P);
			Normalize_by_B_Box_2d<_T>(uv, n, 1,pNorm_uv, NULL, K_uv_Inv);
		}
		Gen_H_Coeff_z_0(pNorm_P, pNorm_uv, n, A);
	}else
		Gen_H_Coeff_z_0(P, uv, n, A);

	//******************接着解齐次矛盾方程组 Ax =0*****************/
	int iResult;
	if (bUse_SVD)//用SVD方法算
		Solve_Homo_Linear_SVD(A, n*2, 9, H, &iResult);
	else
	{
		_T AtA[9 * 9];
		Transpose_Multiply(A, n * 2, 9, AtA, 0);
		iResult = bInverse_Power<_T>(AtA, 9, (_T*)NULL, H, (_T)1e-8);
	}
	if (!iResult)
		goto END;
	//******************接着解齐次矛盾方程组 Ax =0*****************/
	
	if (bNormalize)
	{//此时需要恢复Scale
		if(Norm_Method==Normalize_Method::Dev)
		{
			//Normalize by Dev
			_T K_uv_Inv_1[3 * 3] = { K_uv_Inv[0],0,K_uv_Inv[2],
								0,K_uv_Inv[1],K_uv_Inv[3],
								0,0,1 };
			_T K_P_1[3 * 3] = { K_P[0],0,K_P[2],
								0,K_P[1],K_P[3],
								0,0,1 };
			Matrix_Multiply_3x3<_T>(K_uv_Inv_1, H, H);
			Matrix_Multiply_3x3<_T>(H, K_P_1, H);
		}else if(Norm_Method == Normalize_Method::Hartley)
		{
			//Hartley 方法
			_T K_uv_Inv_1[3 * 3] = { K_uv_Inv[0],0,K_uv_Inv[1],
								0,K_uv_Inv[0],K_uv_Inv[2],
								0,0,1 };
			_T K_P_1[3 * 3] = { K_P[0],0,K_P[1],
								0,K_P[0],K_P[2],
								0,0,1 };
			Matrix_Multiply_3x3<_T>(K_uv_Inv_1, H, H);
			Matrix_Multiply_3x3<_T>(H, K_P_1, H);
		}else if (Norm_Method == Normalize_Method::Bounding_Box)
		{
			//Bounding Box 方法
			_T K_uv_Inv_1[3 * 3] = { K_uv_Inv[0],0,K_uv_Inv[1],
								0,K_uv_Inv[0],K_uv_Inv[2],
								0,0,1 };
			_T K_P_1[3 * 3] = { K_P[0],0,K_P[1],
								0,K_P[0],K_P[2],
								0,0,1 };
			Matrix_Multiply_3x3<_T>(K_uv_Inv_1, H, H);
			Matrix_Multiply_3x3<_T>(H, K_P_1, H);
		}
	}

	bRet = 1;
END:
	Free(A);
	return bRet;
}
template<typename _T>_T Test_H(_T P[][3], _T uv[][2], int n, _T H[3 * 3])
{
	_T fError = 0;
	for (int i = 0; i < n; i++)
	{
		_T uv_1[3];
		Matrix_Multiply(H, 3, 3, P[i], 1, uv_1);
		uv_1[0] /= uv_1[2], uv_1[1] /= uv_1[2];
		fError += (uv_1[0] - uv[i][0]) * (uv_1[0] - uv[i][0]) +
			(uv_1[1] - uv[i][1]) * (uv_1[1] - uv[i][1]);
	}
	fError = (_T)sqrt(fError) / n;
	return fError;
}

template float Test_H_2d(float P[][2], float uv[][2], int n, float H[3 * 3]);
template double Test_H_2d(double P[][2], double uv[][2], int n, double H[3 * 3]);
template<typename _T>_T Test_H_2d(_T P[][2], _T uv[][2], int n, _T H[3 * 3])
{//测试 uv = 1/s * H * P 的误差
	_T(*pP1)[3] = (_T(*)[3])pMalloc(n * 3 * sizeof(_T));
	for (int i = 0; i < n; i++)
		pP1[i][0] = P[i][0], pP1[i][1] = P[i][1], pP1[i][2] = 1;
	_T fError = Test_H(pP1, uv, n, H);
	Free(pP1);
	return fError;
	//return 0;
}

//**********************一组旋转转换************************/
template<typename _T>void Rotation_Vector_3_2_Matrix(_T V[3], _T R[3 * 3])
{
	_T V1[4];
	Rotation_Vector_3_2_4(V, V1);
	Rotation_Vector_4_2_Matrix(V1, R);
}

#define Scale_Matrix_1(A1,fMax) \
{ \
	A1[0][0] *= fMax; \
	A1[0][1] *= fMax; \
	A1[0][2] *= fMax; \
	A1[1][0] *= fMax; \
	A1[1][1] *= fMax; \
	A1[1][2] *= fMax; \
	A1[2][0] *= fMax; \
	A1[2][1] *= fMax; \
	A1[2][2] *= fMax; \
}

template int bIs_Hat(float A[], int n);
template int bIs_Hat(double A[], int n);
template<typename _T>int bIs_Hat(_T A[],int n)
{//验证一个矩阵是否为反对称
	_T* pA1 = (_T*)pMalloc(n * n * sizeof(_T));
	int bRet = 0;
	Matrix_Transpose(A, n, n, pA1);
	Vector_Multiply<_T>(pA1, n * n, -1, pA1);
	_T eps = (_T)1e-10;
	for (int i = 0; i < n * n; i++)
	{
		_T fDelta = pA1[i] - A[i];
		if (Abs(fDelta) < eps)
		{
			Free(pA1);
			return 0;
		}
	}
	Free(pA1);
	return 1;
}


template void Rotation_Vector_4_2_Matrix(float V[4], float R[3 * 3]);
template void Rotation_Vector_4_2_Matrix(double V[4], double R[3 * 3]);
template<typename _T>void Rotation_Vector_4_2_Matrix(_T V[4], _T R[3 * 3])
{//旋转向量到旋转矩阵，试一把看看准不准, 3阶已经一摸一样
//向量的前三个分量为标准化旋转轴，最后一个分量为旋转角
//严格右手系，xyz，x向左，y向远,z向上
	//第一条公式，R = exp(V^)
	////第二条公式，罗德里格斯公式，避开算无限项，此处用第二种
	//_T fCos_Theta, fSin_Theta;
	//_T fTheta = V[3];
	//_T V_1[3];
	////fTheta需要来个求模？
	//fCos_Theta = (_T)cos(fTheta);
	//fSin_Theta = (_T)sin(fTheta);
	//_T nnt[3][3], I[3][3] = { {1,0,0},{0,1,0},{0,0,1} };
	//_T Skew_Sym[3][3];

	////保险起见，规格化一下
	//Normalize(V, 3, V_1);

	//Matrix_Multiply(V_1, 3, 1, V_1, 3, (_T*)nnt);
	//Scale_Matrix_1(I, fCos_Theta);
	//Scale_Matrix_1(nnt, (1 - fCos_Theta));
	//Hat(V_1, (_T*)Skew_Sym);
	//Scale_Matrix_1(Skew_Sym, fSin_Theta);

	//Matrix_Add((_T*)I, (_T*)nnt, 3, (_T*)R);
	//Matrix_Add((_T*)R, (_T*)Skew_Sym, 3, (_T*)R);

	//第三条公式，直接计算
	//展开快一倍，初步具有商用价值
	_T x = V[0], y = V[1], z = V[2],
		cos_theta = (_T)cos(V[3]),
		sin_theta = (_T)sin(V[3]),
		one_minus_cos_theta = (_T)(1 - cos(V[3]));

	//先算对角线
	R[0] = cos_theta + x * x * one_minus_cos_theta;
	R[4] = cos_theta + y * y * one_minus_cos_theta;
	R[8] = cos_theta + z * z * one_minus_cos_theta;

	{	//R01，R101
		_T z_sin_theta = z * sin_theta;
		_T xy_1_cost_theta = x*y * one_minus_cos_theta;
		R[1] = xy_1_cost_theta - z_sin_theta;
		R[3] = xy_1_cost_theta + z_sin_theta;
	}

	{//R02, R20
		_T y_sin_theta = y * sin_theta;
		_T xz_1_cost_theta = x*z * one_minus_cos_theta;
		R[2] = xz_1_cost_theta + y_sin_theta;
		R[6] = xz_1_cost_theta - y_sin_theta;
	}

	{//R12, R21
		_T x_sin_theta = x * sin_theta;
		_T yz_1_cost_theta = y*z * one_minus_cos_theta;
		R[5] = yz_1_cost_theta - x_sin_theta;
		R[7] = yz_1_cost_theta + x_sin_theta;
	}

	return;
}
template void Rotation_Vector_4_2_3(float V[], float V1[]);
template void Rotation_Vector_4_2_3(double V[], double V1[]);
template<typename _T>void Rotation_Vector_4_2_3(_T V[], _T V1[])
{//将4维旋转向量转换维3维旋转向量
	//for (int i = 0; i < 3; i++)
		//V1[i] = V[i] * V[3];
	//展开更快
	V1[0] = V[0] * V[3];
	V1[1] = V[1] * V[3];
	V1[2] = V[2] * V[3];
}
template<typename _T>void Rotation_Vector_3_2_4(_T V[], _T V1[])
{//将3维旋转向量转化为4维旋转向量
	_T fMod = (_T)sqrt(V[0] * V[0] + V[1] * V[1] + V[2] * V[2]);
	V1[0] = V[0] / fMod, V1[1] = V[1] / fMod, V1[2] = V[2] / fMod;
	V1[3] = fMod;
}

template void Rotation_Matrix_2_Vector(float R[3 * 3], float V[4]);
template void Rotation_Matrix_2_Vector(double R[3 * 3], double V[4]);
template<typename _T>void Rotation_Matrix_2_Vector(_T R[3 * 3], _T V[4])
{//从旋转矩阵到旋转向量就是解 Rn=n，其中n就是待求的转轴，显然特征值为1， 求解特征方程
//这个函数有问题
	//由于特征值=1， 代入(A-rI)x=0, 求得x便是特征向量。而r=1,所以解(A-I)x=0即可
	//_T I[3][3] = { 1,0,0,0,1,0,0,0,1 };
	printf("tended to be obsolete\n");

	_T R_1[3][3];// B[3] = { 0 }, V_1[3 * 3]
	_T B[3] = { 0 };
	int iResult;
	if (R == V)
	{
		printf("R connot be V\n");
		return;
	}

	memcpy(R_1, R, 9 * sizeof(_T));
	R_1[0][0] -= 1.f, R_1[1][1] -= 1.f, R_1[2][2] -= 1.f;

	//所以应该用svd分解
	SVD_Info oSVD;
	SVD_Alloc<_T>(3, 3, &oSVD);
	svd_3((_T*)R_1, oSVD, &iResult);

	//Vt的最后一行就是解
	memcpy(V, &((_T*)oSVD.Vt)[6], 3 * sizeof(_T));
	Free_SVD(&oSVD);

	//再求旋转角度，感觉来个负数才行，待考，具体还得验证
	_T fTr = fGet_Tr(R, 3);
	const _T eps = (_T)1e-5;

	_T fTemp = (fTr - 1.f) / 2.f;
	fTemp = Clip3(-1.f, 1.f, fTemp);
	fTemp =  - (_T)acos(fTemp);
	V[3] = fTemp;

	//还有个问题尚未弄利索，旋转向量的正负号问题。这是因为特征向量可正可符，因为
	//特征向量乘以任意常数依旧是原矩阵的特征向量。故此这个向量的正负号是否影响后续
	//的求解，尚待深化
	return;
}
template void Rotation_Matrix_2_Vector_3(float R[3 * 3], float V[3]);
template void Rotation_Matrix_2_Vector_3(double R[3 * 3], double V[3]);
template<typename _T>void Rotation_Matrix_2_Vector_3(_T R[3 * 3], _T V[3])
{//看来书上的解方程不靠谱，因为特征向量可以有无数个，即使模长一样也有两个方向相仿
	_T fCos_Theta = (_T)((R[0] + R[4] + R[8] - 1) * 0.5);
	_T fTheta = (_T)acos(fCos_Theta);
	_T fFactor = (_T)(0.5 * fTheta / sin(fTheta));

	V[0] = (R[2 * 3 + 1] - R[1 * 3 + 2]) * fFactor;
	V[1] = (R[0 * 3 + 2] - R[2 * 3 + 0]) * fFactor;
	V[2] = (R[1 * 3 + 0] - R[0 * 3 + 1]) * fFactor;

	return;
}
template<typename _T>void Rotation_Vector_2_Quaternion(_T V[4], _T Q[4])
{//旋转向量转换为四元组，本来按照定义，一个旋转向量由一个标准化向量作为旋转轴与一个旋转角度构成
	_T fSin_Theta_Div_2;
	_T V_1[3];
	_T fTheta = V[3];
	Normalize(V, 3, V_1);
	Q[0] = cos(fTheta / 2.f);
	fSin_Theta_Div_2 = sin(fTheta / 2.f);
	Q[1] = V_1[0] * fSin_Theta_Div_2;
	Q[2] = V_1[1] * fSin_Theta_Div_2;
	Q[3] = V_1[2] * fSin_Theta_Div_2;
	return;
}
template<typename _T>void Rotation_Matrix_2_Quaternion(_T R[], _T Q[])
{//旋转矩阵到四元组。此处不完全实现，由于缺乏直接算法，故此靠一个旋转向量作为中间商倒腾过去
	_T V[4];
	Rotation_Matrix_2_Vector(R, V);
	Rotation_Vector_2_Quaternion(V, Q);
	return;
}
void Quaternion_Add(float Q_1[], float Q_2[], float Q_3[])
{
	for (int i = 0; i < 4; i++)
		Q_3[i] = Q_1[i] + Q_2[i];
	return;
}
void Quaternion_Minus(float Q_1[], float Q_2[], float Q_3[])
{
	for (int i = 0; i < 4; i++)
		Q_3[i] = Q_1[i] - Q_2[i];
	return;
}
void Quaternion_Conj(float Q_1[], float Q_2[])
{//简单求个共轭
	Q_2[0] = Q_1[0];
	Q_2[1] = -Q_1[1];
	Q_2[2] = -Q_1[2];
	Q_2[3] = -Q_1[3];
}
void Quaternion_Multiply(float Q_1[], float Q_2[], float Q_3[])
{//乘法既不是点积也不是外积，有其定义
	Q_3[0] = Q_1[0] * Q_2[0] - Q_1[1] * Q_2[1] - Q_1[2] * Q_2[2] - Q_1[3] * Q_2[3];
	Q_3[1] = Q_1[0] * Q_2[1] + Q_1[1] * Q_2[0] + Q_1[2] * Q_2[3] - Q_1[3] * Q_2[2];
	Q_3[2] = Q_1[0] * Q_2[2] - Q_1[1] * Q_2[3] + Q_1[2] * Q_2[0] + Q_1[3] * Q_2[1];
	Q_3[3] = Q_1[0] * Q_2[3] + Q_1[1] * Q_2[2] - Q_1[2] * Q_2[1] + Q_1[3] * Q_2[0];
	return;
}
void Quaternion_Inv(float Q_1[], float Q_2[])
{//对四元数求逆
	float fMod = fGet_Mod(Q_1, 4);
	int i;
	Quaternion_Conj(Q_1, Q_2);
	fMod *= fMod;
	for (i = 0; i < 4; i++)
		Q_2[i] /= fMod;
	return;
}
template<typename _T>void Quaternion_2_Rotation_Matrix(_T Q[4], _T R[])
{//四元数转换为旋转矩阵， R= vv' + s^2*I + 2sv^ + (v^)^2
	_T fValue, M_2[3][3], M_1[3][3] = { {1,0,0},{0,1,0},{0,0,1} };	//临时矩阵	
	//先算个vv' 直接放R即可
	Matrix_Multiply(&Q[1], 3, 1, &Q[1], 3, R);

	//再算 s^2*I
	fValue = Q[0] * Q[0];
	Scale_Matrix_1(M_1, fValue);
	Matrix_Add(R, (_T*)M_1, 3, R);

	//再算2sv
	Hat(&Q[1], (_T*)M_1);
	fValue = 2.f * Q[0];
	memcpy(M_2, M_1, 3 * 3 * sizeof(_T));
	Scale_Matrix_1(M_2, fValue);
	Matrix_Add(R, (_T*)M_2, 3, R);

	//再算(v^) ^ 2，上面已经搞定了M_1= v^
	Matrix_Multiply((_T*)M_1, 3, 3, (_T*)M_1, 3, (_T*)M_2);
	Matrix_Add(R, (_T*)M_2, 3, R);

	//Disp((float*)R, 3, 3);
	return;
}
template<typename _T>void Quaternion_2_Rotation_Vector(_T Q[4], _T V[4])
{//由四元数转换为旋转向量,注意了，四元数必须是标准化，即|v|=1，否则旋转角度就不对
	_T fSin_Theta_Div_2;
	V[3] = 2.f * acos(Q[0]);
	fSin_Theta_Div_2 = sin(V[3] / 2.f);
	V[0] = Q[1] / fSin_Theta_Div_2;
	V[1] = Q[2] / fSin_Theta_Div_2;
	V[2] = Q[3] / fSin_Theta_Div_2;
	return;
}
//**********************一组旋转转换************************/

//************************一组投影函数************************/
template void Get_Distort_Coeff(float x, float y, float D[5], float dc[2 * 5]);
template void Get_Distort_Coeff(double x, double y, double D[5], double dc[2 * 5]);
template<typename _T>void Get_Distort_Coeff(_T x, _T y, _T D[5], _T dc[2 * 5])
{//得畸变参数系数矩阵 2x5
	//畸变系数矩阵, Distort Coefficients
	//x1 * r² 	x1 * r^4		x1 * r^6		2*x1 * y1		(r² + 2x1²) 
	//y1 * r²	y1 * r^4 		y1 * r^6		(r² + 2y²) 		2 * x1 * y1
	//Pd = Md * D		2x5 * 5x1 => 2x1 
	_T r2 = x * x + y * y, r4 = r2 * r2, r6 = r2 * r4;
	_T dc_1[] = { x * r2, x * r4,	x * r6,	2 * x * y,	r2 + 2 * x * x,
		y * r2,	y * r4, y * r6,	r2 + 2 * y * y, 2 * x * y };
	memcpy(dc, dc_1, 2 * 5 * sizeof(_T));
	//dc[0] = x * r2, dc[1] = x * r4, dc[2] = x * r6, dc[3] = 2 * x * y, dc[4] = r2 + 2 * x * x;
	//dc[5] = y * r2, dc[6] = y * r4, dc[7] = y * r6, dc[8] = r2 + 2 * y * y, dc[9] = 2 * x * y;
	return;
}

template void Get_Distort_Value(double x, double y, double D[5], double d[2]);
template<typename _T>void Get_Distort_Value(_T x, _T y, _T D[5], _T d[2])
{//对于归一化平面上坐标(x,y) 及给定得畸变参数，求畸变具体数值
//注意，这个只有理论意义，实际很慢
	_T dc[2 * 5];
	Get_Distort_Coeff(x, y, D, dc);
	Matrix_Multiply(dc, 2, 5, D, 1, d);
	return;
}

template void Get_dTP_dKsi(float Pt[3], float Deriv[4 * 6]);
template void Get_dTP_dKsi(double Pt[3], double Deriv[4 * 6]);
template<typename _T>void Get_dTP_dKsi(_T Pt[3], _T Deriv[4 * 6])
{//对于Pt = TP, 求T 上的扰动对Pt的影响
// = dTP/dKsi = dPt/dKsi
	_T P1_M[3 * 3];

	//dTP/dksi = dP'/dksi= I -P'^
	Hat(Pt, P1_M);

	//I
	Deriv[0 * 6 + 0] = 1; Deriv[0 * 6 + 1] = 0; Deriv[0 * 6 + 2] = 0;
	Deriv[1 * 6 + 0] = 0; Deriv[1 * 6 + 1] = 1; Deriv[1 * 6 + 2] = 0;
	Deriv[2 * 6 + 0] = 0; Deriv[2 * 6 + 1] = 0; Deriv[2 * 6 + 2] = 1;
	//-P'^
	Deriv[0 * 6 + 3] = -P1_M[0]; Deriv[0 * 6 + 4] = -P1_M[1]; Deriv[0 * 6 + 5] = -P1_M[2];
	Deriv[1 * 6 + 3] = -P1_M[3]; Deriv[1 * 6 + 4] = -P1_M[4]; Deriv[1 * 6 + 5] = -P1_M[5];
	Deriv[2 * 6 + 3] = -P1_M[6]; Deriv[2 * 6 + 4] = -P1_M[7]; Deriv[2 * 6 + 5] = -P1_M[8];

	//以下只具有理论意义，一般用不上，所以注掉
	//memset(&Deriv[3 * 6], 0, 6 * sizeof(_T));
	return;
}

template void Get_dTP_dKsi(float T[3 * 4], float P[2], float Deriv[4 * 6]);
template void Get_dTP_dKsi(double T[3 * 4], double P[2], double Deriv[4 * 6]);
template<typename _T>void Get_dTP_dKsi(_T T[3*4], _T P[2], _T Deriv[4 * 6])
{//这个是符合函数，先求P' = TP, 再求 dTP/dKsi 然后求导
	_T P1[4] = { P[0],P[1],P[2],1 }, TP[4];
	Matrix_Multiply(T, 3, 4, P1, 1, TP);
	Get_dTP_dKsi(TP, Deriv);
	return;
}
template void Get_dRP_dphi(double Rp[3], double dRP_dphi[3 * 3]);
template<typename _T>void Get_dRP_dphi(_T Rp[3], _T dRP_dphi[3 * 3])
{//Rp对李代数phi 的求导。李导数上的左绕顶对Rp 位置的影响
	//= -(Rp)^
	_T hat[3 * 3];
	Hat(Rp, hat);
	Vector_Multiply<_T>(hat, 3 * 3, -1, dRP_dphi);
}

template void Get_dPn_dP(double Pt[3], double R[3 * 3], double dPn_dP[2 * 3]);
template<typename _T>void Get_dPn_dP(_T Pt[3], _T R[3 * 3], _T dPn_dP[2 * 3])
{//空间点（世界坐标下）P 的微小扰动对归一化平面上投影的影响
	//其中 Pt 为TP 变换后的相机坐标，R 为位姿T的旋转
	//dPn/dPt = dPn/dPt * dPt/dP
	/*_T dPn_dPt[2 * 3];
	Get_dPn_dPt(Pt, R, dPn_dPt);
	Matrix_Multiply(dPn_dPt, 2, 3, R, 3, dPn_dP);*/

	//最终展开最快
	_T fRecip = 1 / Pt[2];
	dPn_dP[0] = fRecip * (R[0 * 3 + 0] - Pt[0] * R[2 * 3 + 0]);
	dPn_dP[1] = fRecip * (R[0 * 3 + 1] - Pt[0] * R[2 * 3 + 1]);
	dPn_dP[2] = fRecip * (R[0 * 3 + 2] - Pt[0] * fRecip * R[2 * 3 + 2]);

	dPn_dP[3] = fRecip * (R[1 * 3 + 0] - Pt[1] * R[2 * 3 + 0]);
	dPn_dP[4] = fRecip * (R[1 * 3 + 1] - Pt[1] * R[2 * 3 + 1]);
	dPn_dP[5] = fRecip * (R[1 * 3 + 2] - Pt[1] * fRecip * R[2 * 3 + 2]);
}

template void Get_dPn_dPt(double Pt[3], double R[3 * 3], double dPn_dPt[2 * 3]);
template<typename _T>void Get_dPn_dPt(_T Pt[3], _T R[3*3], _T dPn_dPt[2*3])
{// Pn = Pt/z, 求 dPn/dPt，这就是	投影雅可比矩阵
//真正意义是相机坐标下的点Pt 微小扰动对归一化平面上投影的影响
	dPn_dPt[0] = dPn_dPt[4] = 1 / Pt[2];
	dPn_dPt[2] = -Pt[0] * dPn_dPt[0] * dPn_dPt[0];
	dPn_dPt[5] = -Pt[1] * dPn_dPt[0] * dPn_dPt[0];
	dPn_dPt[1] = dPn_dPt[3] = 0;
}

template void Get_dPd_dPn(double xn, double yn, double k1, double k2, double k3, double p1, double p2, double dPd_dPn[4]);
template<typename _T>void Get_dPd_dPn(_T xn, _T yn, _T k1, _T k2, _T k3, _T p1, _T p2, _T dPd_dPn[4])
{//(x,y) 为归一化平面未畸变前的坐标，畸变是求导中的难点，特别容易出错
	_T r2 = xn * xn + yn * yn, r4 = r2 * r2, r6 = r2 * r4;
	//1 + k1 * r² + k2 * r^4 + k3 * r^6
	_T d = 1 + k1 * r2 + k2 * r4 + k3 * r6;
	//dd / dr2 = k1 + 2 * k2 * r2 + 3 * K3 * r2 ^ 2
	_T dd_dr2 = k1 + 2 * k2 * r2 + 3 * k3 * r2 * r2;

	//d + 2x^2  * dd/dr2 + 2*p1*y + 6*p2*x
	dPd_dPn[0] = d + 2 * xn * xn * dd_dr2 + 2 * p1 * yn + 6 * p2 * xn;
	//2*x*y * dd/dr2 + 2*p1*x + 2*p2*y
	dPd_dPn[2] = dPd_dPn[1] = 2 * xn * yn * dd_dr2 + 2 * p1 * xn + 2 * p2 * yn;
	//d + 2y^2 * dyd/dr2 + 6 * p1 * y + 2*p2*x 
	dPd_dPn[3] = d + 2 * yn * yn * dd_dr2 + 6 * p1 * yn + 2 * p2 * xn;
}

template void Get_PnP_Deriv(float P[4], float uv_Ref[2], float T[4 * 4], float K[3 * 3], float D[5], float dE_dK[], float dE_dKsi[], float dE_dD[], float dE_dP[], float E[2]);
template void Get_PnP_Deriv(double P[4], double uv_Ref[2], double T[4 * 4], double K[3 * 3], double D[5], double dE_dK[], double dE_dKsi[], double dE_dD[], double dE_dP[], double E[2]);
template<typename _T>void Get_PnP_Deriv(_T P[4], _T uv_Ref[2],_T T[4 * 4], _T K[3 * 3], _T D[5],
	_T dE_dK[],_T dE_dKsi[], _T dE_dD[], _T dE_dP[],_T E[2])
{//搞个满血版的求导，只作为基准程序对数据用
	//第一步，求TP
	_T Pt[3];
	Matrix_Multiply(T, 3, 4, P, 1, Pt);
	//Disp(Pt, 1, 3, "Pt");

	//投影到归一化平面
	_T x = Pt[0] / Pt[2], y = Pt[1] / Pt[2];
	_T Pn[2] = { x,y };

	//根据畸变参数与归一化坐标算出具体畸变得距离
	_T d[2], dc[2 * 5];
	Get_Distort_Coeff(x, y, D, dc);
	Matrix_Multiply(dc, 2, 5, D, 1, d);

	_T Pd[2];
	Vector_Add(Pn, d, 2, Pd);

	//再将归一化平面投影到像素平面
	_T UV[2];
	UV[0] = Pd[0] * K[0] + K[2];
	UV[1] = Pd[1] * K[4] + K[5];
	/*Disp(K, 3, 3, "K");
	Disp(D, 1, 5, "D");
	Disp(T, 4, 4, "T");
	Disp(P, 1, 3, "P");*/
	//Disp(UV, 1, 2, "uv");

	_T E1[2];
	E1[0] = uv_Ref[0] - UV[0];
	E1[1] = uv_Ref[1] - UV[1];
	//Disp(uv_Ref, 1, 2, "uv'");
	//Disp(E1, 1, 2, "E");

	//轮到求导
	_T dLoss_dE[2] = { -1,-1 };

	_T dLoss_dUV[2] = { -E1[0],-E1[1] };	//dLoss/dE = (-eu,	-ev)

	_T dE_dK1[2 * 4] = { -Pd[0],	0, -1,	0,
						0,		-Pd[1],		0, -1 };
	//Disp(dE_dK1, 2, 4, "dE/dK");

	//对畸变后的坐标进行求导
	_T fx = K[0], fy = K[4];
	_T dE_dPd[2 * 2] = { -fx,0,
						0,-fy };
	//Disp(dE_dPd, 2, 2, "dE/dPd");

	_T dE_dD1[2 * 5];	//= dE/dPd * dt
	Matrix_Multiply(dE_dPd, 2, 2, dc, 5, dE_dD1);
	//Disp(dE_dPd, 2, 2, "dE/dPd");
	//Disp(dE_dD1, 2, 5, "dE/dD");

	_T dPd_dPn[4];
	Get_dPd_dPn(x, y, D[0], D[1], D[2], D[3], D[4], dPd_dPn);
	//Disp(dPd_dPn, 2, 2, "dPd/dPn");

	_T dE_dPn[2*2];	//dE/dPn = dE/dPd * dPd/dPn 	2x2 * 2x2 => 2x2
	Matrix_Multiply(dE_dPd, 2, 2, dPd_dPn, 2, dE_dPn);
	//Disp(dE_dPn, 2, 2, "dE/dPn");

	_T Ptz_Recip = 1 / Pt[2], Ptz_Recip_Sqr = Ptz_Recip * Ptz_Recip;
	//投影雅可比矩阵
	_T dPn_dPt[2 * 3] = {	Ptz_Recip,	0,			-Pt[0] * Ptz_Recip_Sqr,
							0,			Ptz_Recip ,	-Pt[1] * Ptz_Recip_Sqr };
	//Disp(dPn_dPt, 2, 3, "dPn/dPt");

	_T dE_dPt[2 * 3];
	Matrix_Multiply(dE_dPn, 2, 2, dPn_dPt, 3, dE_dPt);
	//Disp(dE_dPt, 2, 3, "dE/dPt");

	_T dTP_dKsi[4 * 6];
	Get_dTP_dKsi(Pt, dTP_dKsi);
	//Disp(dTP_dKsi, 3, 6, "dTP/dKsi");

	_T dE_dKsi1[2 * 6];
	Matrix_Multiply(dE_dPt, 2, 3, dTP_dKsi, 6, dE_dKsi1);
	//Disp(dE_dKsi1, 2, 6, "dE/dKsi");

	_T dE_dP1[2 * 3];	//= dE/dPt * R
	_T R[3 * 3];
	Get_R_t<_T>(T, R, NULL);
	Matrix_Multiply(dE_dPt, 2, 3, R, 3, dE_dP1);
	//Disp(dE_dP1, 2, 3, "dE/dP");

	if (dE_dK)
		memcpy(dE_dK, dE_dK1, 2 * 4 * sizeof(_T));
	if (dE_dKsi)
		memcpy(dE_dKsi, dE_dKsi1, 2 * 6 * sizeof(_T));
	if (dE_dD)
		memcpy(dE_dD, dE_dD1, 2 * 5 * sizeof(_T));
	if (dE_dP)
		memcpy(dE_dP, dE_dP1, 2 * 3 * sizeof(_T));
	if (E)
		memcpy(E, E1, 2 * sizeof(_T));

	return;
}

template<typename _T>void Get_Pd_by_Pn(_T Pn[2], _T D[5], _T Pd[])
{//根据归一化平面位置，与畸变参数，求畸变后位置
	/*_T Delta[2];
	Get_Distort_Value<_T>(Pn[0], Pn[1], D, Delta);
	Vector_Add(Pn, Delta, 2, Pd);*/

	//1 + k1 * r2 + k2 * r^4 + k3 * r^6)
	_T x = Pn[0], y = Pn[1];
	_T r2 = x * x + y * y, r4 = r2 * r2, r6 = r2 * r4;
	_T d = 1 + D[0] * r2 + D[1] * r4 + D[2] * r6;
	_T xy = x * y;
	Pd[0] = x * d + 2 * D[3] * xy + D[4] * (r2 + 2 * x * x);
	Pd[1] = y * d + D[3] * (r2 + 2 * y * y) + 2 * D[4] * xy;
	return;
}

template<typename _T>static void Get_Pn_by_Pd_JtJ_JtE(_T J[2 * 2], _T E[2],
	_T JtJ[5], _T JtE[2])
{//精简版展开得 J'J J'E
	//J'J
	JtJ[0] = J[0] * J[0] + J[1] * J[1];
	JtJ[1] = JtJ[2] = J[0] * J[1] + J[1] * J[3];
	JtJ[3] = J[1] * J[1] + J[3] * J[3];

	//J'E
	JtE[0] = J[0] * E[0] + J[2] * E[1];
	JtE[1] = J[1] * E[0] + J[3] * E[1];
	return;
}

template int Get_Pn_by_Pd(double Pd[2], double Distort[5], double Pn[2]);
template<typename _T>int  Get_Pn_by_Pd(_T Pd[2], _T Distort[5], _T Pn[2])
{//用畸变值xd 反推归一化平面上的真实值
	_T fError, fPre_Error = MAX_FLOAT, fError_Delta, eps = 1e-15;
	const int iMax_Iter = 30;
	int bRet = 0, iIter, iResult;

	//设置求解参数初值
	_T x[2] = { Pd[0],Pd[1] }, Pd_1[2];
	for (iIter = 0; iIter < iMax_Iter; iIter++)
	{//迭代体
		//注意：以下将Pd 作为Pn 的初值
		Get_Pd_by_Pn<_T>(x, Distort, Pd_1);
		_T  E[2] = { Pd[0] - Pd_1[0],Pd[1] - Pd_1[1] };
		fError = fGet_Sqr_Sum(E, 2);
		fError_Delta = fPre_Error - fError;
		if (fError < eps)
			break;
		if (fError_Delta < 0)
			goto END;   //发散了
		fPre_Error = fError;

		_T J[2 * 2];
		Get_dPd_dPn(x[0], x[1], Distort[0], Distort[1], Distort[2], Distort[3], Distort[4], J);

		_T JtJ[4], JtE[2];    //精简版上三角系数矩阵
		Get_Pn_by_Pd_JtJ_JtE(J, E, JtJ, JtE);  //通过J,E生成 Ax = b

		//解方程 Jx = E
		_T Delta_x[2];
		if (!(iResult = Solve_Linear_2x2(JtJ, JtE, Delta_x)))
			goto END;

		x[0] += Delta_x[0], x[1] += Delta_x[1];
	}
	if (iIter >= iMax_Iter)
		goto END;

	Pn[0] = x[0], Pn[1] = x[1];
	bRet = 1;
END:
	return bRet;
}

template void Get_uv_Ref(float P[4], float T[4 * 4], float K[3 * 3], float D[5], float uv[2]);
template void Get_uv_Ref(double P[4], double T[4 * 4], double K[3 * 3], double D[5], double uv[2]);
template<typename _T>void Get_uv_Ref(_T P[4], _T T[4 * 4], _T K[3 * 3], _T D[5], _T uv[2])
{//搞一个满血版的从空间点开始，经过相机投影，畸变，投影到像素平面
	//第一步，求TP
	_T Pt[3];
	Matrix_Multiply(T, 3, 4, P, 1, Pt);
	//Disp(Pt, 1, 3, "Pt");
	//投影到归一化平面
	_T x = Pt[0] / Pt[2], y = Pt[1] / Pt[2];
	_T Pn[2] = { x,y };
	//Disp(Pn, 1, 2, "Pn");
	//根据畸变参数与归一化坐标算出具体畸变得距离
	_T d[2];
	Get_Distort_Value(x, y, D, d);

	_T Pd[2];
	Vector_Add(Pn, d, 2, Pd);
	//Disp(Pd, 1, 2, "Pd");

	//再将归一化平面投影到像素平面
	uv[0] = Pd[0] * K[0] + K[2];
	uv[1] = Pd[1] * K[4] + K[5];
	//printf("u= %f*%f+%f= %f", Pd[0], K[0], K[2], uv[0]);

	return;
}
//************************一组投影函数************************/

//***********************李群李代数**********************************/
template void Get_J_E_uv(float P[4], float uv[2], float T[3 * 4], float K[4], float D[5], float J[2][15], float E[2]);
template void Get_J_E_uv(double P[4], double uv[2], double T[3 * 4], double K[4], double D[5], double J[2][15], double E[2]);
template<typename _T>void Get_J_E_uv(_T P[4], _T uv[2], _T T[3 * 4], _T K[4], _T D[5], _T J[2][15], _T E[2])
{//在像素平面上求雅可比，残差
///给定一个空间点，给定相机位姿T, 内参K,D，一个对应的uv，算个雅可比，误差E
	_T dE_dKsi[2 * 6], dE_dK[2 * 4], dE_dD[2 * 5];
	Get_PnP_Deriv<_T>(P, uv, T, K, D, dE_dK, dE_dKsi, dE_dD, NULL, E);
	//Disp(dE_dK, 2, 4, "dE/dK");
	//Disp(dE_dD, 2, 5, "dE/dD");
	//Disp(dE_dKsi, 2, 6, "dE/dKsi");
	Copy_Matrix_Partial<_T>(dE_dKsi, 2, 6, (_T*)J, 15, 0, 0);
	Copy_Matrix_Partial<_T>(dE_dK, 2, 4, (_T*)J, 15, 6, 0);
	Copy_Matrix_Partial<_T>(dE_dD, 2, 5, (_T*)J, 15, 10, 0);
	return;
}

template void TP(double T[3 * 4], double P[3], double Pt[3]);
template<typename _T>void TP(_T T[3 * 4], _T P[3], _T Pt[3])
{//Pt: 表示P经过相机位姿T 变换后的相机位置 Pt = TP
	if (P == Pt)
	{//这就不安全了
		_T Pt1[3];
		Pt1[0] = T[0] * P[0] + T[1] * P[1] + T[2] * P[2] + T[3];
		Pt1[1] = T[4] * P[0] + T[5] * P[1] + T[6] * P[2] + T[7];
		Pt1[2] = T[8] * P[0] + T[9] * P[1] + T[10] * P[2] + T[11];
		memcpy(Pt, Pt1, 3 * sizeof(_T));
	}else
	{
		Pt[0] = T[0] * P[0] + T[1] * P[1] + T[2] * P[2] + T[3];
		Pt[1] = T[4] * P[0] + T[5] * P[1] + T[6] * P[2] + T[7];
		Pt[2] = T[8] * P[0] + T[9] * P[1] + T[10] * P[2] + T[11];
	}
}

template<typename _T>void Get_PnP_Norm_Deriv(_T P[4], _T uv_Ref[2], _T T[4 * 4], _T K[3 * 3], _T D[5],
	_T dE_dK[2 * 4] = NULL, _T dE_dKsi[2 * 6] = NULL, _T dE_dD[2 * 5] = NULL, _T dE_dP[2 * 3] = NULL, _T E[2] = NULL)
{
	_T Pt[3];
	TP(T, P, Pt);

	//投影到归一化平面
	_T x = Pt[0] / Pt[2], y = Pt[1] / Pt[2];
	_T Pn[2] = { x,y };
	//Disp(Pn, 1, 2, "Pn");

	_T d[2], dc[2 * 5];
	Get_Distort_Coeff<_T>(x, y, D, dc);
	Matrix_Multiply(dc, 2, 5, D, 1, d);

	_T Pd[2];
	Vector_Add(Pn, d, 2, Pd);
	//Disp(Pd, 1, 2, "Pd");

	_T E1[2], UVn[2] = { (uv_Ref[0] - K[2]) / K[0], (uv_Ref[1] - K[3]) / K[1] };
	Vector_Minus(UVn, Pd, 2, E1);
	//Disp(E1, 1, 2, "E");

	_T dE_dK1[] = { -(uv_Ref[0] - K[2]) / (K[0] * K[0]),0, -1 / K[0],0,
					0,-(uv_Ref[1] - K[3]) / (K[1] * K[1]),0,-1 / K[1] };
	//Disp(dE_dK1, 2, 4, "dE/dK");

	////对归一化平面坐标进行求导
	//_T r2 = x * x + y * y, r4 = r2 * r2, r6 = r4 * r2;
	//_T dn = 1 + D[0] * r2 + D[1] * r4 + D[2] * r6;
	//_T dPd_dPn[4] = { dn + 2 * D[3] * Pn[1] + 4 * D[4] * Pn[0] ,	2 * D[3] * Pn[1],
	//			2 * D[4] * Pn[1],	dn + 4.f * D[3] * Pn[1] + 2 * D[4] * Pn[0] };

	_T dPd_dPn[4];
	Get_dPd_dPn(x, y, D[0], D[1], D[2], D[3], D[4], dPd_dPn);

	_T dE_dPd[2 * 2] = { -1,0,
				0,-1 };
	//Disp(dE_dPd, 2, 2, "dE/dPd");

	_T dE_dD1[2 * 5];	//= dE/dPd * dt
	Matrix_Multiply(dE_dPd, 2, 2, dc, 5, dE_dD1);
	//Disp(dE_dD1, 2, 5, "dE/dD");
	_T dE_dPn[2 * 2];
	Matrix_Multiply(dE_dPd, 2, 2, dPd_dPn, 2, dE_dPn);

	_T Ptz_Recip = 1 / Pt[2], Ptz_Recip_Sqr = Ptz_Recip * Ptz_Recip;
	_T dPn_dPt[2 * 3] = { Ptz_Recip,0,-Pt[0] * Ptz_Recip_Sqr,
					0,Ptz_Recip ,-Pt[1] * Ptz_Recip_Sqr };

	_T dE_dPt[2 * 3];
	Matrix_Multiply(dE_dPn, 2, 2, dPn_dPt, 3, dE_dPt);
	//Disp(dE_dPt, 2, 3, "dE/dPt");

	_T dTP_dKsi[4 * 6];
	Get_dTP_dKsi(Pt, dTP_dKsi);
	Disp(dTP_dKsi, 3, 6, "dTP/dKsi");

	_T dE_dKsi1[2 * 6];
	Matrix_Multiply(dE_dPt, 2, 3, dTP_dKsi, 6, dE_dKsi1);
	//Disp(dE_dKsi1, 2, 6, "dE/dKsi");

	_T dE_dP1[2 * 3];	//= dE/dPt * R
	_T R[3 * 3];
	Get_R_t<_T>(T, R, NULL);
	Matrix_Multiply(dE_dPt, 2, 3, R, 3, dE_dP1);
	//Disp(dE_dP1, 2, 3, "dE/dP");

	if (dE_dK)
		memcpy(dE_dK, dE_dK1, 2 * 4 * sizeof(_T));
	if (dE_dKsi)
		memcpy(dE_dKsi, dE_dKsi1, 2 * 6 * sizeof(_T));
	if (dE_dD)
		memcpy(dE_dD, dE_dD1, 2 * 5 * sizeof(_T));
	if (dE_dP)
		memcpy(dE_dP, dE_dP1, 2 * 3 * sizeof(_T));
	if (E)
		memcpy(E, E1, 2 * sizeof(_T));

	return;
}

template void Get_J_E_Norm(float P[4], float uv[2], float T[3 * 4], float K[4], float D[5], float J[2][15], float E[2]);
template void Get_J_E_Norm(double P[4], double uv[2], double T[3 * 4], double K[4], double D[5], double J[2][15], double E[2]);
template<typename _T>void Get_J_E_Norm(_T P[4], _T uv[2], _T T[3 * 4], _T K[4], _T D[5], _T J[2][15], _T E[2])
{//在归一化平面上求雅可比，残差
///给定一个空间点，给定相机位姿T, 内参K,D，一个对应的uv，算个雅可比，误差E
	_T dE_dKsi[2 * 6], dE_dK[2 * 4], dE_dD[2 * 5];
	Get_PnP_Norm_Deriv<_T>(P, uv, T, K, D, dE_dK, dE_dKsi, dE_dD, NULL, E);
	Copy_Matrix_Partial<_T>(dE_dKsi, 2, 6, (_T*)J, 15, 0, 0);
	Copy_Matrix_Partial<_T>(dE_dK, 2, 4, (_T*)J, 15, 6, 0);
	Copy_Matrix_Partial<_T>(dE_dD, 2, 5, (_T*)J, 15, 10, 0);
	return;
}

template void se3_2_SE3(float Ksi[6], float T[]);
template void se3_2_SE3(double Ksi[6], double T[]);
template<typename _T>void se3_2_SE3(_T Ksi[6], _T T[])
{//T是6维 se3向量Ksi对应的4x4矩阵, Ksi前rho后phi
//注意：感觉这个转坏是错误的，问题在J上
//转换完毕后，T完全与图形学的三维转换矩阵一致
//总结， SE3中的旋转矩阵，se3中的六维向量，都是先旋转后位移
	//首先求R
	_T R[3][3];
	_T Rotation_Vector[4];

	Rotation_Vector_3_2_4(&Ksi[3], Rotation_Vector);
	Rotation_Vector_4_2_Matrix(Rotation_Vector, (_T*)R);
	//Disp((_T*)R, 3, 3,"R");

	_T J[3][3], J_Rho[3];

	//问题可能就在下面
	//显然，J与ρ无关，只从φ推导出来
	Get_Jl_4(Rotation_Vector, (_T*)J);
	Matrix_Multiply((_T*)J, 3, 3, Ksi, 1, J_Rho);

	//然后将 R,J_Rho, 0', 1组合成T
	T[0] = R[0][0], T[1] = R[0][1], T[2] = R[0][2], T[3] = J_Rho[0];
	T[4] = R[1][0], T[5] = R[1][1], T[6] = R[1][2], T[7] = J_Rho[1];
	T[8] = R[2][0], T[9] = R[2][1], T[10] = R[2][2], T[11] = J_Rho[2];
	T[12] = T[13] = T[14] = 0, T[15] = 1;

	return;
}

template void Gen_Pose_By_V3_t(float R[], float t[], float T[]);
template void Gen_Pose_By_V3_t(double R[], double t[], double T[]);
template<typename _T>void Gen_Pose_By_V3_t(_T V3[], _T t[], _T T[])
{
	_T R[3 * 3];
	Rotation_Vector_3_2_Matrix(V3, R);
	Gen_Pose_By_R_t(R, t, T);
}
template void Gen_Pose_By_R_t(float R[], float t[], float T[], int b4x4);
template void Gen_Pose_By_R_t(double R[], double t[], double T[], int b4x4);
template<typename _T>void Gen_Pose_By_R_t(_T R[], _T t[], _T T[], int b4x4)
{//安全函数，目标可以等于源
//用旋转坐标与位移坐标构成一个4x4 齐次变换矩阵，此处由旋转与平移构成
//这个矩阵的物理意义应该是先旋转后平移
	_T T1[4 * 4];
	if (R)
	{
		T1[0] = R[0], T1[1] = R[1], T1[2] = R[2];
		T1[4] = R[3], T1[5] = R[4], T1[6] = R[5];
		T1[8] = R[6], T1[9] = R[7], T1[10] = R[8];
	}
	else
	{//此处代考
		memset(T1, 0, 4 * 4 * sizeof(_T));
		T1[0] = T[1 * 4 + 1] = T[2 * 4 + 2] = 0;
	}
	T1[15] = 1;
	if (t)
	{
		T1[3] = t[0];
		T1[7] = t[1];
		T1[11] = t[2];
	}
	else
		T1[3] = T1[7] = T1[11] = 0;
	if(b4x4)
	{
		T1[12] = T1[13] = T1[14] = 0;
		memcpy(T, T1, 4 * 4 * sizeof(_T));
	}else
		memcpy(T, T1, 3 * 4 * sizeof(_T));

	return;
}

template void Get_Inv_T(float T[3 * 4], float T_Inv[3 * 4], int b4x4);
template void Get_Inv_T(double T[3 * 4], double T_Inv[3 * 4], int b4x4);
template<typename _T>void Get_Inv_T(_T T[3 * 4], _T T_Inv[3 * 4], int b4x4)
{//利用T = Rt 的物理性质，驯熟求得T 的逆矩阵
	_T R[3 * 3] = { T[0],T[4],T[8],		//利用正交矩阵性质，直接转置了事
				T[1],T[5],T[9],
				T[2],T[6],T[10] };
	_T t[3] = { -T[3],-T[7],-T[11] };
	Matrix_Multiply_3x1(R, t, t);
	Gen_Pose_By_R_t(R, t, T_Inv, b4x4);
}

template void T_2_R9_t(float T[3 * 4], float R[3 * 3], float t[3]);
template void T_2_R9_t(double T[3 * 4], double R[3 * 3], double t[3]);
template<typename _T>void T_2_R9_t(_T T[3 * 4], _T R[3 * 3], _T t[3])
{
	Copy_Matrix_Partial(T, 4, 4, R, 3, 0, 0);
	t[0] = T[3], t[1] = T[7], t[2] = T[11];
}
template void Vee(float M[], float V[3]);
template void Vee(double M[], double V[3]);
template<typename _T>void Vee(_T M[], _T V[3])
{//反对称矩阵到向量
	V[0] = M[7];
	V[1] = M[2];
	V[2] = M[3];
	return;
}

template void Hat(float V[], float M[]);
template void Hat(double V[], double M[]);
template<typename _T>void Hat(_T V[], _T M[])
{//根据给定的向量构造反对称矩阵，改回与书中一致
	if (V == M)
	{
		_T M1[3 * 3];
		M1[0] = M1[4] = M1[8] = 0;
		M1[1] = -V[2], M1[3] = V[2];
		M1[2] = V[1], M1[6] = -V[1];
		M1[5] = -V[0], M1[7] = V[0];
		memcpy(M, M1, 3 * 3 * sizeof(_T));
	}
	else
	{
		M[0] = M[4] = M[8] = 0;
		M[1] = -V[2], M[3] = V[2];
		M[2] = V[1], M[6] = -V[1];
		M[5] = -V[0], M[7] = V[0];
	}
	return;
}
template void Get_R_t(float T[4 * 4], float R[3 * 3], float t[3]);
template void Get_R_t(double T[4 * 4], double R[3 * 3], double t[3]);
template<typename _T>void Get_R_t(_T T[4 * 4], _T R[3 * 3], _T t[3])
{//从4x4 齐次矩阵中抽取R，t
	if (R)
	{
		R[0] = T[0], R[1] = T[1], R[2] = T[2];
		R[3] = T[4], R[4] = T[5], R[5] = T[6];
		R[6] = T[8], R[7] = T[9], R[8] = T[10];
	}
	if (t)
		t[0] = T[3], t[1] = T[7], t[2] = T[11];
	return;
}

template void Get_Jr_4(float phi[4], float J[]);
template void Get_Jr_4(double phi[4], double J[]);
template<typename _T>void Get_Jr_4(_T phi[4], _T J[])
{
	_T v[4];
	Vector_Multiply<_T>(phi, 4, -1, v);
	Get_Jl_4(v, J);
}

template void Get_Jl_4(float phi[4],  float J[]);
template void Get_Jl_4(double phi[4],  double J[]);
template<typename _T>void Get_Jl_4(_T phi[4], _T J[])
{//给定的旋转向量，求出J矩阵。旋转向量为4元组
	//再求位移t,先求出个J，a为旋转向量的转轴
	_T fValue, Temp_1[3][3], I[3][3] = { {1,0,0},{0,1,0},{0,0,1} };
	int i;
	memset(J, 0, 3 * 3 * sizeof(_T));

	//先求J第一部分 (sin(theta)/theta)*I
	if (phi[3] != 0)
		fValue = (_T)sin(phi[3]) / phi[3];
	else
		fValue = 0;
	for (i = 0; i < 9; i++)
		((_T*)J)[i] = fValue * ((_T*)I)[i];

	//再求J第二部分 (1-sin(theta)/theta) * axa'
	fValue = 1.f - fValue;
	Matrix_Multiply(phi, 3, 1, phi, 3, (_T*)Temp_1);
	//Disp((float*)Temp_1, 3, 3);
	for (i = 0; i < 9; i++)
		((_T*)J)[i] += fValue * ((_T*)Temp_1)[i];

	//再求第三部分 (1-cos(theta))/theta * a^
	if (phi[3] != 0)
		fValue = (_T)(1 - cos(phi[3])) / phi[3];
	else
		fValue = 0;
	Hat(phi, (_T*)Temp_1);
	for (i = 0; i < 9; i++)
		((_T*)J)[i] += fValue * ((_T*)Temp_1)[i];
}

template void Get_Jl_3(float phi[3], float J[]);
template void Get_Jl_3(double phi[3], double J[]);
template<typename _T>void Get_Jl_3(_T phi[3], _T J[])
{//李代数上的扰动引起李群的变化
	_T v[4];
	Rotation_Vector_3_2_4(phi, v);
	Get_Jl_4(v, J);
	Matrix_Transpose(J, 3, 3, J);
}

template void Get_Jr_3(float phi[3], float J[]);
template void Get_Jr_3(double phi[3], double J[]);
template<typename _T>void Get_Jr_3(_T phi[3], _T J[])
{//李代数上的扰动引起李群的变化
	_T v[4];
	Rotation_Vector_3_2_4(phi, v);
	Vector_Multiply<_T>(v, 4, -1, v);
	Get_Jl_4(v, J);
	Matrix_Transpose(J, 3, 3, J);
}

template void Exp_V3(float V[], float B[], float eps);
template void Exp_V3(double V[], double B[], double eps);
template<typename _T>void Exp_V3(_T V[3], _T B[383], _T eps)
{
	_T hat[3 * 3];
	Hat(V, hat);
	Exp_M(hat, 3, B, eps);
}

//***********************李群李代数**********************************/

//******************************画出位姿*****************************************************/
template void Draw_Camera(Point_Cloud<float>* poPC, float T[4 * 4], int R, int G, int B);
template void Draw_Camera(Point_Cloud<double>* poPC, double T[4 * 4], int R, int G, int B);
template<typename _T>void Draw_Camera(Point_Cloud<_T>* poPC, _T T[4 * 4], int R, int G, int B)
{//看看能否画出个简陋的相机位姿

	//先对T求逆
	_T T_Inv[4 * 4];
	//int iResult;
	Get_Inv_T(T, T_Inv);

	_T View_Point[4] = { T_Inv[3],T_Inv[7],T_Inv[11],1 };    //此处应该算是视点
	_T Norm_Plane[4][4] = { {-1,1,-1,1 },
		{1,1,-1,1},
		{1,-1,-1,1},
		{-1,-1,-1,1} };    //归一化平面上的4个点
	_T Norm_Center[4] = { 0,0,-1,1 };      //归一化平面上的中心
	//_T Temp[4];
	int i;
	//剩下的一律从原地搬过去
	//Draw_Sphere(poPC, View_Point[0], View_Point[1], View_Point[2],(_T)0.2f,40,255,0,0);
	Draw_Point(poPC, View_Point[0], View_Point[1], View_Point[2], R, G, B);

	Matrix_Multiply(T_Inv, 4, 4, Norm_Center, 1, Norm_Center);
	Draw_Line(poPC, View_Point[0], View_Point[1], View_Point[2], Norm_Center[0], Norm_Center[1], Norm_Center[2], 50, R, G, B);

	for (i = 0; i < 4; i++)
	{
		Matrix_Multiply(T_Inv, 4, 4, Norm_Plane[i], 1, Norm_Plane[i]);
		Draw_Line(poPC, View_Point[0], View_Point[1], View_Point[2], Norm_Plane[i][0], Norm_Plane[i][1], Norm_Plane[i][2], 50, R, G, B);
	}
	for (i = 0; i < 4; i++)
		Draw_Line(poPC, Norm_Plane[i][0], Norm_Plane[i][1], Norm_Plane[i][2], Norm_Plane[(i + 1) & 3][0], Norm_Plane[(i + 1) & 3][1], Norm_Plane[(i + 1) & 3][2], 50, R, G, B);

	Draw_Line(poPC, Norm_Plane[0][0], Norm_Plane[0][1], Norm_Plane[0][2], Norm_Plane[2][0], Norm_Plane[2][1], Norm_Plane[2][2], 50, R, G, B);
	Draw_Line(poPC, Norm_Plane[1][0], Norm_Plane[1][1], Norm_Plane[1][2], Norm_Plane[3][0], Norm_Plane[3][1], Norm_Plane[3][2], 50, R, G, B);
	return;
}
//******************************画出位姿*****************************************************/

template<typename _T>void Get_H_Block(_T J[2][15], _T E[2], _T Camera[6 * 6], _T Cam_Corner[6 * 9], _T Corner[9 * 9], _T JtE[6 + 9])
{//单独拎出来搞块
	//Camera	属于位姿块，6x6	左上角
	//Corner	属于内参块，9x9	右下角
	//Cam_Corner 属于位姿-内参块,右上，左下角，堆成，所以只寸一个
	//JtE:		也要累加进去
	int i, j, iDest_Pos;

	////造一些数据来看看
	//for (int i = 0; i < 2; i++)
	//	for (int j = 0; j < 15; j++)
	//		J[i][j] = (i + 1) + j;

	//_T Temp[15 * 15] = {};
	//Transpose_Multiply<_T>((_T*)J, 2, 15, Temp, 0);

	//Disp((_T*)J, 2, 15, "J");
	//Disp(Temp, 15, 15, "Temp");

	for (i = 0; i < 6; i++)
	{
		iDest_Pos = i * 6 + i;
		for (j = i; j < 6; j++, iDest_Pos++)
		{
			//iDest_Pos = i * 6 + j;
			Camera[iDest_Pos] += J[0][i] * J[0][j] + J[1][i] * J[1][j];
			//Corner[iDest_Pos] += J[0][i + 6] * J[0][j + 6] + J[1][i + 6] * J[1][j + 6];
		}
	}
	//Disp(Camera, 6, 6, "Camera");
	for (i = 0; i < 9; i++)
	{
		iDest_Pos = i * 9 + i;
		for (j = i; j < 9; j++, iDest_Pos++)
			Corner[iDest_Pos] += J[0][i + 6] * J[0][j + 6] + J[1][i + 6] * J[1][j + 6];
	}
	//Disp(Corner, 9, 9, "Corner");
	for (i = 0; i < 6; i++)
		for (j = 0; j < 9; j++)
			Cam_Corner[i * 9 + j] += J[0][i] * J[0][j + 6] + J[1][i] * J[1][j + 6];

	for (int i = 0; i < 15; i++)
		JtE[i] += J[0][i] * E[0] + J[1][i] * E[1];

	//Disp(Cam_Corner, 6, 9, "Cam_Corner");
	return;
}

template void Get_H_JtE(float T[][3 * 4], float K[4], float D[5], float P[][2], Point_2D<float> uv[], int iObservation_Count, int iOrder, float H[], float JtE[]);
template void Get_H_JtE(double T[][3 * 4], double K[4], double D[5], double P[][2], Point_2D<double> uv[], int iObservation_Count, int iOrder, double H[], double JtE[]);
template<typename _T>void Get_H_JtE(_T T[][3 * 4], _T K[4], _T D[5], _T P[][2], Point_2D<_T> uv[],	int iObservation_Count, int iOrder, _T H[], _T JtE[])
{//给定所有的位姿，K,D，生成一个H = J'J, 预计J'E
	_T Camera[6 * 6] = {}, Cam_KD[6 * 9] = {}, K_D[9 * 9] = {},
		JtE_1[6 + 9] = {};
	int iPre_Cam_Index = uv[0].m_iCamera_Index,
		iCount_Minus_1 = iObservation_Count - 1;
	memset(JtE, 0, iOrder * sizeof(_T));
	memset(H, 0, iOrder * iOrder * sizeof(_T));

	_T K9[3 * 3];
	K4_2_K9(K, K9);
	for (int i = 0; i < iObservation_Count;)
	{
		//每一点都能算出一个雅可比，E
		Point_2D<_T> oUV = uv[i];
		if (oUV.m_iCamera_Index != iPre_Cam_Index)
		{
			//补下三角
			for (int y = 1; y < 6; y++)
			{
				for (int x = 0; x < y; x++)
				{
					int iDest_Pos = y * 6 + x,
						iSource_Pos = x * 6 + y;
					Camera[iDest_Pos] = Camera[iSource_Pos];
				}
			}for (int y = 1; y < 9; y++)
			{
				for (int x = 0; x < y; x++)
				{
					int iDest_Pos = y * 9 + x,
						iSource_Pos = x * 9 + y;
					K_D[iDest_Pos] = K_D[iSource_Pos];
				}
			}

			int iCam_Index = iPre_Cam_Index != -1 ? iPre_Cam_Index : oUV.m_iCamera_Index;

			//将块拷到sigma_H
			Copy_Matrix_Partial(Camera, 6, 6, H, iOrder, iCam_Index * 6, iCam_Index * 6);
			Copy_Matrix_Partial(Cam_KD, 6, 9, H, iOrder, iOrder - 9, iCam_Index * 6);

			memcpy(JtE + iCam_Index * 6, JtE_1, 6 * sizeof(_T));
			Vector_Add(JtE + iOrder - 9, JtE_1 + 6, 9, JtE + iOrder - 9);

			_T* pDest_1 = &H[(iOrder - 9) * iOrder + iCam_Index * 6];
			for (int y = 0; y < 9; y++, pDest_1 += iOrder)
				for (int x = 0; x < 6; x++)
					pDest_1[x] = Cam_KD[x * 9 + y];

			if (iPre_Cam_Index == -1)
				break;
			memset(Camera, 0, 6 * 6 * sizeof(_T));
			memset(Cam_KD, 0, 6 * 9 * sizeof(_T));
			memset(JtE_1, 0, 15 * sizeof(_T));
			iPre_Cam_Index = oUV.m_iCamera_Index;
			continue;
		}

		_T J[2][6 + 9], E[2];
		_T* P1 = P[oUV.m_iPoint_Index];
		_T P2[4] = { P1[0],P1[1],0,1 };

		Get_J_E_uv(P2, oUV.m_Pos, T[oUV.m_iCamera_Index], K9, D, J, E);

		//然后用J 算一个 H = J'J
		Get_H_Block<_T>(J, E, Camera, Cam_KD, K_D, JtE_1);

		if (i == iCount_Minus_1)
			iPre_Cam_Index = -1;
		else
			i++;
	}
	//Disp_Fillness(H, iWidth_H, iWidth_H, "H");
	Copy_Matrix_Partial(K_D, 9, 9, H, iOrder, iOrder - 9, iOrder - 9);

	////验算
	//memset(H, 0, 15 * 15 * sizeof(_T));
	//memset(JtE, 0, 15 * sizeof(_T));
	//_T fTotal = 0;
	//for (int i = 0; i < iObservation_Count;i++)
	//{
	//	Point_2D<_T> oUV = uv[i];
	//	_T J[2][6 + 9], E[2], JtJ[15 * 15];
	//	_T* P1 = P[oUV.m_iPoint_Index];
	//	_T P2[3] = { P1[0],P1[1],0 };
	//	Get_J_E_Norm(P2, oUV.m_Pos, T[oUV.m_iCamera_Index], K, D, J, E);
	//	//Disp((_T*)J, 2, 15, "J");
	//	//Disp((_T*)E, 2, 1, "E");

	//	Transpose_Multiply((_T*)J, 2, 15, JtJ, 0);
	//	//printf("%f\n", JtJ[0]);
	//	Vector_Add(H, JtJ, 15 * 15, H);
	//	//fTotal += JtJ[0];
	//	//printf("fTotal:%f H[0]:%f\n", fTotal, H[0]);

	//	At_x_B((_T*)J, 2, 15, E, 1, JtE_1);
	//	printf("%f\n", JtE_1[0]);
	//	fTotal += JtE_1[0];
	//	Vector_Add(JtE, JtE_1, 15, JtE);
	//}
	////printf("%f\n", fGet_Cond_Num(H, 15));
	////Disp(H, 15, 15, "H");
	////Disp(JtE, 15, 1, "JtE");
	return;
}

//*********************************求E矩阵**************************/
template<typename _T>static void Gen_E_xEy_row(_T x[2], _T y[2], _T row[9])
{//xEy = 0 方程的系数矩阵的一行   咬死 x2' * E * x1 = 0
	row[0] = x[0] * y[0];
	row[1] = x[1] * y[0];
	row[2] = y[0];
	row[3] = x[0] * y[1];
	row[4] = x[1] * y[1];
	row[5] = y[1];
	row[6] = x[0];
	row[7] = x[1];
	row[8] = 1;

	//以下是颠倒顺序，x1' * E * x2 = 0
	/*row[0] = x[0] * y[0];
	row[1] = x[0] * y[1];
	row[2] = x[0];
	row[3] = x[1] * y[0];
	row[4] = x[1] * y[1];
	row[5] = x[1];
	row[6] = y[0];
	row[7] = y[1];
	row[8] = 1;*/
	return;
}

template<typename _T>static _T Get_E_Poly_Coeff(_T E_Base[4][3 * 3], int a, int b, int c)
{//构造一个矩阵用于算行列式
	_T E[3 * 3];
	memcpy(E, E_Base[a], 3 * sizeof(_T));
	memcpy(&E[3], &E_Base[b][1 * 3], 3 * sizeof(_T));
	memcpy(&E[6], &E_Base[c][2 * 3], 3 * sizeof(_T));
	//Disp(E, 3, 3, "E");
	return fGet_Determinant<_T>(E, 3);
}

template<typename _T>static _T Contribute(_T E_Base[4][3 * 3], _T fab[4][4][3][3], int r, int c, int a, int b, int d)
{//r,c 为Mrc，最终矩阵的位置
	_T fCon, fPart_1 = 0, fPart_2 = 0;
	for (int l = 0; l < 3; l++)
	{
		fPart_1 += fab[a][b][r][l] * E_Base[d][l * 3 + c];
		fPart_2 += fab[a][b][l][l];
	}

	fPart_1 *= 2;
	fPart_2 *= E_Base[d][r * 3 + c];
	fCon = fPart_1 - fPart_2;
	return fCon;
}

template<typename _T>static void Get_E_Det_Coeff(_T E_Base[4][3 * 3], _T Coeff[20])
{
	//先构造行列式约束
	_T Map[4][4][4];
	for (int a = 0; a < 4; a++)
	{
		_T* ra = E_Base[a];
		for (int b = 0; b < 4; b++)
		{
			_T* rb = &E_Base[b][1 * 3];
			//显然，以下是余子式
			_T r01_12_02_11 = ra[1] * rb[2] - ra[2] * rb[1],
				r00_12_02_00 = ra[0] * rb[2] - ra[2] * rb[0],
				r00_11_01_10 = ra[0] * rb[1] - ra[1] * rb[0];

			for (int c = 0; c < 4; c++)
			{
				_T* rc = &E_Base[c][2 * 3];
				//显然，以下求出了行列式
				Map[a][b][c] = rc[0] * r01_12_02_11 - rc[1] * r00_12_02_00 +
					rc[2] * r00_11_01_10;
				//Map[a][b][c] = Get_E_Poly_Coeff(E_Base, a, b, c);
			}
		}
	}

	//x开头项
	//x^3
	Coeff[0] = Map[0][0][0];
	//x^2*y
	Coeff[1] = Map[0][0][1] + Map[0][1][0] + Map[1][0][0];
	//x^2 * z
	Coeff[2] = Map[0][0][2] + Map[0][2][0] + Map[2][0][0];
	//x^2
	Coeff[3] = Map[0][0][3] + Map[0][3][0] + Map[3][0][0];
	//x * y^2
	Coeff[4] = Map[0][1][1] + Map[1][0][1] + Map[1][1][0];
	//xyz
	Coeff[5] = Map[0][1][2] + Map[0][2][1] + Map[1][0][2] + Map[1][2][0] + Map[2][0][1] + Map[2][1][0];
	//x*y
	Coeff[6] = Map[0][1][3] + Map[0][3][1] + Map[1][0][3] + Map[1][3][0] + Map[3][0][1] + Map[3][1][0];
	//x*z^2
	Coeff[7] = Map[0][2][2] + Map[2][0][2] + Map[2][2][0];
	//x*z
	Coeff[8] = Map[0][2][3] + Map[0][3][2] + Map[2][0][3] + Map[2][3][0] + Map[3][0][2] + Map[3][2][0];
	//x
	Coeff[9] = Map[0][3][3] + Map[3][0][3] + Map[3][3][0];

	//y开头项
	//y^3
	Coeff[10] = Map[1][1][1];
	//y^2 * z
	Coeff[11] = Map[1][1][2] + Map[1][2][1] + Map[2][1][1];
	//y^2
	Coeff[12] = Map[1][1][3] + Map[1][3][1] + Map[3][1][1];
	//y*z^2
	Coeff[13] = Map[1][2][2] + Map[2][1][2] + Map[2][2][1];
	//y*z
	Coeff[14] = Map[1][2][3] + Map[1][3][2] + Map[2][1][3] + Map[2][3][1] + Map[3][1][2] + Map[3][2][1];
	//y
	Coeff[15] = Map[1][3][3] + Map[3][1][3] + Map[3][3][1];
	//Disp(&Coeff[10], 6, 1, "y");

	//z^3
	Coeff[16] = Map[2][2][2];
	//z^2
	Coeff[17] = Map[2][2][3] + Map[2][3][2] + Map[3][2][2];
	//z
	Coeff[18] = Map[2][3][3] + Map[3][2][3] + Map[3][3][2];

	//常数项目
	Coeff[19] = Map[3][3][3];
	//Disp(Coeff, 20, 1, "Coeff");
	return;
}

template<typename _T>static void Get_E_Trace_Coeff(_T E_Base[4][3 * 3], _T Coeff[9][20])
{//算迹约束方程系数
	_T fab[4][4][3][3];
	for (int a = 0; a < 4; a++)
		for (int b = 0; b < 4; b++)
			for (int r1 = 0; r1 < 3; r1++)
			{
				_T* pr1 = &E_Base[a][r1 * 3];
				for (int r2 = 0; r2 < 3; r2++)
				{//抽出a矩阵r1 行
					_T* pr2 = &E_Base[b][r2 * 3];
					fab[a][b][r1][r2] = fDot(pr1, pr2, 3);
				}
			}

	const char Exp_2_Index[20][3] = {
				3,0,0,  //x^3
				2,1,0,  //x^2 * y
				2,0,1,  //x^2 * z 
				2,0,0,  //x^2   
				1,2,0,  //x * y^2     
				1,1,1,  //xyz
				1,1,0,  //x*y     
				1,0,2,  //x*z^2
				1,0,1,  //x*z     
				1,0,0,  //x

				0,3,0,  //y^3   
				0,2,1,  //y^2 * z
				0,2,0,  //y^2
				0,1,2,  //y*z^2
				0,1,1,  //y*z
				0,1,0,  //y

				0,0,3,  //z^3
				0,0,2,  //z^2
				0,0,1,  //z

				0,0,0   //常数项
	};

	for (int r = 0; r < 3; r++)
	{
		for (int c = 0; c < 3; c++)
		{
			_T Exp_444[4][4][4] = {};    //装累加项
			for (int a = 0; a < 4; a++)
				for (int b = 0; b < 4; b++)
					for (int d = 0; d < 4; d++)
					{
						unsigned char Exp_4[4] = {};
						Exp_4[a]++;
						Exp_4[b]++;
						Exp_4[d]++;
						_T fValue = Contribute(E_Base, fab, r, c, a, b, d);
						Exp_444[Exp_4[0]][Exp_4[1]][Exp_4[2]] += fValue;
					}

			_T* pCoeff = Coeff[r * 3 + c];
			for (int i = 0; i < 20; i++)
			{
				const char* pMap = Exp_2_Index[i];
				pCoeff[i] = Exp_444[pMap[0]][pMap[1]][pMap[2]];
			}
		}
	}
	return;
}

template<typename _T>static void Get_E_Coeff_10x20(_T E_Base[4][3 * 3], _T Coeff[10 * 20])
{
	Get_E_Det_Coeff(E_Base, Coeff);
	Get_E_Trace_Coeff(E_Base, (_T(*)[20]) & Coeff[20]);
}

template<typename _T>static void Re_Arrange_E_Coeff(_T Coeff[10 * 20])
{//当行列式约束与迹约束系数矩阵求出后，本来按照x,y,z 顺序排，
   //要改成按次数由高到低排
	_T Temp[10 * 20];
	unsigned char Map[20] = { 0,1,2,4,5,7,10,11,13,16,3,6,8,12,14,17,9,15,18,19 };
	for (int y = 0; y < 10; y++)
	{
		_T* pOrg = &Coeff[y * 20];
		_T* pTemp = &Temp[y * 20];
		for (int x = 0; x < 20; x++)
			pTemp[x] = pOrg[Map[x]];
	}
	memcpy(Coeff, Temp, 10 * 20 * sizeof(_T));
	return;
}

template<typename _T>static void Gen_E_Az(_T Coeff[10 * 20], _T Az[10 * 10])
{//由B 构成早Az 矩阵
	//x3 * Xbase =	x1^2 * x3				Xhigh [2]		-B[2]
	Vector_Multiply<_T>(&Coeff[2 * 20 + 10], 10, -1, &Az[0]);
	//x1*x2 *x3				Xhigh [4]		-B[4]
	Vector_Multiply<_T>(&Coeff[4 * 20 + 10], 10, -1, &Az[1 * 10]);
	//x1*x3^2				Xhigh [5]		-B5]
	Vector_Multiply<_T>(&Coeff[5 * 20 + 10], 10, -1, &Az[2 * 10]);
	//x2^2 * x3				Xhigh [7]		-B[7]
	Vector_Multiply<_T>(&Coeff[7 * 20 + 10], 10, -1, &Az[3 * 10]);
	//x2*x3^2				Xhigh [8]		-B[8]
	Vector_Multiply<_T>(&Coeff[8 * 20 + 10], 10, -1, &Az[4 * 10]);
	//X3 ^ 3				Xhigh[9] - B[9]
	Vector_Multiply<_T>(&Coeff[9 * 20 + 10], 10, -1, &Az[5 * 10]);
	//0 0 1 0 0 0 0 0 0 0 
	_T Aux[10 * 4] = { 0,0,1,0,0,0,0,0,0,0,
						0,0,0,0,1,0,0,0,0,0,
						0,0,0,0,0,1,0,0,0,0,
						0,0,0,0,0,0,0,0,1,0 };
	memcpy(&Az[6 * 10], Aux, 10 * 4 * sizeof(_T));

	return;
}

template<typename _T>static int Get_E_Solution(_T Az[], _T X[10][3], int* piCount)
{//对Az 求解特征值，从z 求出x,y，可能由多组,后续再验根
	Complex<_T> x1[10];
	int iCount, iResult;
	iResult = QR_Decompose_Complex<_T>(Az, 10, x1, &iCount);

	//以下求解z 对应的x,y, 求解路径是先求解z对应的特征向量，
	//此时一个z 对应多个特征向量，每个向量再求出对应的x,y,
	//本来可能会出现一个z对应多个x,y。但是5点法的几何意义注定
	//一个z只能对应一个x,y，古翠不用复杂化
	for (int i = 0; i < iCount; i++)
	{
		//求关于z 的特征向量
		_T V[10];   // , AzV[10];
		int iVector_Count;
		Get_Eigen_Vector_by_Value_Real(Az, 10, x1[i].real, V, &iVector_Count);
		if (iVector_Count > 1)
		{
			printf("Invalid eigen vector count\n");
			return 0;
		}
		if (iVector_Count == 0)
		{
			printf("0 eigen vector\n");
			continue;
		}

		//Matrix_Multiply(Az, 10, 10, V, 1, AzV);
		//Disp(AzV, 10, 1, "AzV");
		//V 对应的向量：
		//      x1^2
		//      x1*x2
		//      x1*x3
		//      x2^2
		//      x2*x3
		//      x3^2
		//      x1
		//      x2
		//      x3
		//      1

		//注意，以下相当于 恢复scale
		//V = s * (x^2, x*y, x*z, y^2, y*z, z^2, x, y, z, 1)'
		//(x^2, x*y, x*z, y^2, y*z, z^2, x, y, z, 1)'/V[9];

		X[i][0] = V[6] / V[9];
		X[i][1] = V[7] / V[9];
		X[i][2] = V[8] / V[9];
	}

	//Disp((_T*)X, iCount, 3,"xyz");
	if (*piCount)
		*piCount = iCount;
	return iResult;
}

template<typename _T>static void E_2_R_t(_T E[], _T R[2][3 * 3], _T t[3])
{//注意，此处求得的t可以直接乘以两相机之间的距离，即可得到真实的位移
	_T U[3 * 3], Vt[3 * 3], S[3 * 3];
	//Disp(E, 3, 3, "E");
	int iResult = SVD_Decompose<_T>(E, 3, S, U, Vt);
	//W =   0	-1	0
	//      1	0	0
	//      0	0	1
	_T W[] = { 0,-1,	0,      //就是 (0,0,PI/2) 对应得旋转矩阵
				1,	0,	0,
				0,	0,	1 };

	_T Temp[3 * 3];

	//生成 R = U * W * V'
	Matrix_Multiply_3x3<_T>(U, W, Temp);
	Matrix_Multiply_3x3<_T>(Temp, Vt, R[0]);

	//t = U*W**S*U'
	memset(S, 0, 9 * sizeof(_T));

	//注意：此处等于将E 归一化。本身从 x1' * E * x2 = 0 可以判断出
	// E 具有尺度不变性。任你将 E放大多少倍，经过 S = {1, 1, 0} 后
	// 全部统一为唯一一个数
	S[0] = S[4] = 1;
	_T S_Ut[3 * 3];
	iResult = A_x_Bt(S, 3, 3, U, 3, S_Ut);
	Matrix_Multiply_3x3(Temp, S_Ut, Temp);
	Vee(Temp, t);
	//Vector_Multiply<_T>(t[0], 3, -1, t[1]);     //t1 = -t0，恒等

	//生成 R = U * W * V'
	W[1] = 1, W[3] = -1;
	Matrix_Multiply_3x3<_T>(U, W, Temp);
	Matrix_Multiply_3x3<_T>(Temp, Vt, R[1]);
	Matrix_Multiply_3x3(Temp, S_Ut, Temp);
	Vee(Temp, t);

	//最后更新E
	memset(S, 0, 9 * sizeof(_T));
	S[0] = S[4] = 1;
	//Disp(S, 3, 3, "S");

	Matrix_Multiply_3x3(U, S, E);
	Matrix_Multiply_3x3(E, Vt, E);
	return;
}

template int Elect_E(double Point_A[5][2], double Point_B[5][2], double E[][3 * 3], int iCount, double Best_T[3 * 4], int iSample_Count, int* piBest_E_Index);
template<typename _T>int Elect_E(_T Point_A[5][2], _T Point_B[5][2], _T E[][3 * 3], int iCount, _T Best_T[3 * 4], int iSample_Count, int* piBest_E_Index)
{//选出一个E
	int bFound = 0, iBest_E_Index = -1;
	_T fError_Sum, fBest_Error = MAX_FLOAT;

	for (int i = 0; i < iCount; i++)
	{
		_T R[2][3 * 3], t0[3], T2[4][4 * 4];
		E_2_R_t(E[i], R, t0);
		_T t1[3] = { -t0[0], -t0[1], -t0[2] };

		Gen_Pose_By_R_t(R[0], t0, T2[0]);
		Gen_Pose_By_R_t(R[0], t1, T2[1]);
		Gen_Pose_By_R_t(R[1], t0, T2[2]);
		Gen_Pose_By_R_t(R[1], t1, T2[3]);
		for (int k = 0; k < 4; k++)
		{
			_T* pT2 = T2[k];
			int bInvalid = 0;
			fError_Sum = 0;

			_T R3[3 * 3], t3[3];    //带个3表示用于三角化
			Get_R_t(pT2, R3, t3);

			for (int j = 0; j < iSample_Count; j++)
			{
				_T P[4];
				Triangulate_Cramer<_T>(Point_A[j], Point_B[j], pT2, P);
				if (P[2] <= 0)
				{
					bInvalid = 1;
					break;
				}
				_T P2[3];
				TP(pT2, P, P2);
				fError_Sum += Test_Triangulate(P, pT2, Point_A[j], Point_B[j]);
				if (P2[2] <= 0)
				{
					bInvalid = 1;
					break;
				}
			}
			if (!bInvalid && fError_Sum < fBest_Error)
			{//成功了
				//printf("E:%d RT:%d\n", i, k);
				if(Best_T)
					memcpy(Best_T, pT2, 3 * 4 * sizeof(_T));
				iBest_E_Index = i;
				fBest_Error = fError_Sum;
				bFound++;
			}
		}
	}

	//注意，bFound 不是所有及格的数量，因为后面的误差比前面大就会掠过不统计
	if (!bFound)
		return 0;
	if (piBest_E_Index)
		*piBest_E_Index = iBest_E_Index;
	return 1;
}

template<typename _T>int Solve_E_5x9_Eq(_T Point_A[5][2], _T Point_B[5][2], _T Root[5 * 9], int* piRoot_Count)
{
	_T A[5 * 9];
	for (int i = 0; i < 5; i++)
		Gen_E_xEy_row(Point_A[i], Point_B[i], &A[i * 9]);

	int iRoot_Count = 0;
	return Solve_Linear_Solution_Construction<_T>(A, 5, 9, NULL, Root, piRoot_Count);
}

template<typename _T>int Estimate_E_5_Point(_T Point_A[5][2], _T Point_B[5][2], _T(**ppE)[3 * 3], int* piCount)
{//硅基出一组E矩阵
	int iResult, bRet = 0;
	int iSize = 10 * 3 * 3 +         //E
		5 * 9 +     //Root, E_Base
		10 * 20 +           //Coeff
		10 * 10 +            //Az
		10 * 3;            //X

	_T* pBuffer = (_T*)pMalloc(iSize * sizeof(_T));
	_T(*pE)[3 * 3] = (_T(*)[3 * 3])pBuffer,
		* pRoot = (_T*)(pE + 10),
		* pCoeff = pRoot + 5 * 9,
		* pAz = pCoeff + 10 * 20,
		(*pX)[3] = (_T(*)[3])(pAz + 10 * 10),
		(*pE_Base)[3 * 3] = NULL;

	//第一步，构造5*9 系数矩阵A
	int iRoot_Count;
	if (!(iResult = Solve_E_5x9_Eq(Point_A, Point_B, pRoot, &iRoot_Count)))
		goto END;
	//Disp(pRoot, 4, 9, "Root");

	//先求出多项式矩阵的系数
	pE_Base = (_T(*)[3 * 3])pRoot;
	Get_E_Coeff_10x20(pE_Base, pCoeff);
	//Disp(pCoeff, 10, 20, "10x20");

	Re_Arrange_E_Coeff(pCoeff);
	//Disp(pCoeff, 10, 20, "10x20");
	if (!(iResult = Elementary_Row_Op_Pivot(pCoeff, 10, 20)))
		goto END;

	//Disp(pCoeff, 10, 20, "10x20");
	//求解行列式与迹约束方程
	Gen_E_Az(pCoeff, pAz);
	//Disp(pAz, 10, 10, "Az");
	if (!(iResult = Get_E_Solution(pAz, pX, &iRoot_Count)))
		goto END;
	//Disp((_T*)pX, iRoot_Count, 3,"xyz");

	//生成一组候选E
	for (int i = 0; i < iRoot_Count; i++)
	{
		_T x = pX[i][0], y = pX[i][1], z = pX[i][2];
		for (int j = 0; j < 3 * 3; j++)
			pE[i][j] = x * pE_Base[0][j] + y * pE_Base[1][j] + z * pE_Base[2][j] + pE_Base[3][j];

		//归一化
		//Normalize(pE[i], 9, pE[i]);
	}

	//缩小已分配内存
	Shrink(pE, iRoot_Count * 3 * 3 * sizeof(_T));
	if (piCount)
		*piCount = iRoot_Count;
	if (ppE)
		*ppE = pE;
	bRet = 1;
END:
	if (!bRet)
	{
		Free(pBuffer);
		pE = NULL;
	}
	return bRet;
}

template int Estimate_E_T_5Point(double Point_A[5][2], double Point_B[5][2], double T[3 * 4], double E[3 * 3]);
template<typename _T>int Estimate_E_T_5Point(_T Point_A[5][2], _T Point_B[5][2], _T T[3 * 4], _T E[3 * 3])
{//求出最后的位姿 T
	//注意，五点法只接受归一化平面上无畸变的点位置，这才符合求解过程
	_T(*pE)[3 * 3] = NULL;
	int iCount, iResult, bRet = 0;
	if (!(iResult = Estimate_E_5_Point<_T>(Point_A, Point_B, &pE, &iCount)))
		goto END;

	//Disp((_T*)pE, iCount, 9, "E");

	//for(int i=0;i<iCount;i++)
		//Test_E(pE[i], Point_A, Point_B,5);

	int iBest_E_Index;
	if(!(iResult = Elect_E<_T>(Point_A, Point_B, pE, iCount, T,5, &iBest_E_Index)))
		goto END;

	if (E)
		memcpy(E, pE[iBest_E_Index], 3 * 3 * sizeof(_T));

	bRet = 1;
END:
	Free(pE);
	return bRet;
}

template int Estimate_E_Rt_5Point(double Point_A[5][2], double Point_B[5][2], double R[3 * 3], double t[3], double E[3 * 3]);
template<typename _T>int Estimate_E_Rt_5Point(_T Point_A[5][2], _T Point_B[5][2], _T R[3 * 3], _T t[3], _T E[3 * 3])
{//求出最后的位姿 R, t
//注意，五点法只接受归一化平面上无畸变的点位置，这才符合求解过程
	_T T[4 * 4];
	int iResult = Estimate_E_T_5Point(Point_A, Point_B, T, E);
	if (!iResult)
		return 0;
	Get_R_t(T, R, t);
	return iResult;
}

template int Estimate_E_T_nPoint(double Point_A[5][2], double Point_B[5][2], int iCount, int* piInlier, double T[3 * 4], double E[3 * 3]);
template<typename _T>int Estimate_E_T_nPoint(_T Point_A[5][2], _T Point_B[5][2], int iCount, int* piInlier, 
	_T T[3 * 4], _T E[3 * 3])
{//用超过8点估计E矩阵，注意了，这是开放算法，不保证所有点对都及格，有些会成为外点，甚至
	//外点比内殿更多

	int bRet = 0, iMax_Inlier = 0, iBest_E_Index = -1;
	_T fError_Sum,* pA = (_T*)pMalloc(iCount * 9 * sizeof(_T));
	if (!pA)
		goto END;
	for(int i=0;i<iCount;i++)
		Gen_E_xEy_row(Point_A[i], Point_B[i], &pA[i * 9]);
	//Disp(pA, iCount, 9, "A");

	int iResult;
	_T x[8 * 9];
	Solve_Linear_Contradictory<_T>(pA, iCount, 9, NULL, x, &iResult);
	if (!iResult)
		return 0;

	//接着需要分解R,t
	_T R[2][3 * 3], t0[3], t1[3], T2[4][4 * 4];
	E_2_R_t(x, R, t0);
	t1[0] = -t0[0], t1[1] = -t0[1], t1[2] = -t0[2];
	Gen_Pose_By_R_t(R[0], t0, T2[0]);
	Gen_Pose_By_R_t(R[0], t1, T2[1]);
	Gen_Pose_By_R_t(R[1], t0, T2[2]);
	Gen_Pose_By_R_t(R[1], t1, T2[3]);
	for (int k = 0; k < 4; k++)
	{
		_T* pT2 = T2[k], P[4];
		int bInvalid = 0, iInlier = 0;
		fError_Sum = 0;

		_T R3[3 * 3], t3[3];    //带个3表示用于三角化
		Get_R_t(pT2, R3, t3);
		
		for (int j = 0; j < iCount; j++)
		{
			Triangulate_Cramer<_T>(Point_A[j], Point_B[j], pT2, P);
			if (P[2] <= 0)
				continue;
			
			_T P2[3];
			TP(pT2, P, P2);
			fError_Sum += Test_Triangulate(P, pT2, Point_A[j], Point_B[j]);
			if (P2[2] <= 0)
				continue;

			iInlier++;
		}
		if (iInlier > iMax_Inlier)
		{
			iMax_Inlier = iInlier;
			iBest_E_Index = k;
		}
	}

	if (E)
		memcpy(E, x, 3 * 3 * sizeof(_T));
	if (T)
		memcpy(T, T2[iBest_E_Index], 3 * 4 * sizeof(_T));
	if (piInlier)
		*piInlier = iMax_Inlier;

	bRet = 1;
END:
	Free(pA);
	return bRet;
}

template int Estimate_E_8Point(double Point_A[8][2], double Point_B[8][2], double T[3 * 4], double E[3 * 3]);
template<typename _T>int Estimate_E_8Point(_T Point_A[8][2], _T Point_B[8][2], _T T[3 * 4], _T E[3 * 3])
{//这个比5点法简单很多
	_T A[8 * 9];
	for (int i = 0; i < 8; i++)
		Gen_E_xEy_row(Point_A[i], Point_B[i], &A[i * 9]);

	//Disp(A, 8, 9, "A");
	int iResult, iSolution_Count;
	_T x[8 * 9];
	iResult = Solve_Linear_Solution_Construction<_T>(A, 8, 9, NULL, x, &iSolution_Count);
	if (!iResult || iSolution_Count != 1)
		return 0;
	//Disp(x, 1, 9, "x");
	int iBest_E_Index;
	if (!(iResult = Elect_E<_T>(Point_A, Point_B, (_T(*)[9])x, 1, T,8, &iBest_E_Index)))
		return 0;
	//Disp(x, 1, 9, "x");
	if (E)
		memcpy(E, x, 3 * 3 * sizeof(_T));
	return 1;
}

template<typename _T>void Estimate_F_Get_F(_T F1[3 * 3], _T F2[3 * 3], _T F[3][3 * 3], int* piCount)
{//从Ax = 0 解除的基构造F
	_T M00[3], M11[3], M01[3];

	//M01 = r1 x s2 + s1 x r2
	Cross_Product(&F1[3], &F2[6], M00);
	Cross_Product(&F2[3], &F1[6], M01);
	Vector_Add(M00, M01, 3, M01);

	//M00 = r1 * r2
	Cross_Product(&F1[3], &F1[6], M00);
	//M11 = s1 * s2
	Cross_Product(&F2[3], &F2[6], M11);

	//a = r0 * M00
	_T a = fDot(F1, M00, 3);
	//b = s0 * M00 + r0 * M01
	_T b = fDot(F2, M00, 3) + fDot(F1, M01, 3);
	//d = r0 * M11 + s0 * M01
	_T d = fDot(F1, M11, 3) + fDot(F2, M01, 3);
	//e = s0 * M11
	_T e = fDot(F2, M11, 3);

	Complex<_T> Root[3];
	int iReal_Count;
	Solve_Cubic(a, b, d, e, Root, &iReal_Count);
	for (int i = 0; i < iReal_Count; i++)
	{
		Vector_Multiply(F1, 3 * 3, Root[i].real, F[i]);
		Vector_Add(F[i], F2, 3 * 3, F[i]);
	}

	//Disp_Complex(Root[i]);
	if (piCount)
		*piCount = iReal_Count;
	return;
}

template int Estimate_F_7Point(double Point_A[7][2], double Point_B[7][2], double F[3][3 * 3], int* piCount);
template<typename _T>int Estimate_F_7Point(_T Point_A[7][2], _T Point_B[7][2], _T F[3][3 * 3], int* piCount)
{//注意，7点法求F矩阵不需要内参，归一化用hartlay
//与E矩阵的求解不一样，此处的求解场合为不知内参。古翠，什么也没有
//由于解三次方程最多有三个根，所以此处最多有3组候选F矩阵
//7点法无法判断哪个F最优，所有判断最优是后续Ransac的事
	//Normalize_Hartley_2d()
	_T Norm_A[7][2], Norm_B[7][2],
		K_Inv_A[4], K_Inv_B[4];
	Normalize_2d<_T>(Point_A, 7, Norm_A, Normalize_Method::Hartley, NULL, K_Inv_A);
	Normalize_2d<_T>(Point_B, 7, Norm_B, Normalize_Method::Hartley, NULL, K_Inv_B);

	//此处沿用x2'*E*x1 = 0 的列式
	union {
		_T A[7 * 9];
		_T x[7 * 9];
	};
	for (int i = 0; i < 7; i++)
		Gen_E_xEy_row(Norm_A[i], Norm_B[i], &A[i * 9]);

	//Disp(A, 7, 9, "A");
	int iResult, iCount;
	iResult = Solve_Linear_Solution_Construction<_T>(A, 7, 9, NULL, x, &iCount);
	if (!iResult || iCount != 2)
		return 0;

	//Disp(x, 2, 3 * 3, "x");
	//_T F_3[3][3 * 3];
	Estimate_F_Get_F(&x[0], &x[9], F, piCount);
	return 1;
}

template int Estimate_F_8Point(double Point_A[8][2], double Point_B[8][2], double F[3 * 3]);
template<typename _T>int Estimate_F_8Point(_T Point_A[8][2], _T Point_B[8][2], _T F[3 * 3])
{//8点法估计F矩阵 x2'*F*x1 =0
	//第一步，Hartlay 归一法
	//注意，归一化不能称为投影，因为投影是降维过程，而归一化唯独
	//不变。此处的变换看上取象内参K投影，实际只是一个相似矩阵
	//本来正路是在归一化平面上做，但是如果归一化平面都有了，还不如
	//用五点法求E矩阵。所以，这个方法解决无内参下求F矩阵的问题
	_T Norm_A[8][2], Norm_B[8][2], T3_A[4], T3_B[4];
	Normalize_2d<_T>(Point_A, 8, Norm_A, Hartley, T3_A);
	Normalize_2d<_T>(Point_B, 8, Norm_B, Hartley, T3_B);

	union {
		_T A[8 * 9];
		_T x[8 * 9];
		_T Temp[3 * 3];
	};
	//求解 x2'*A*x1 = 0
	for (int i = 0; i < 8; i++)
		Gen_E_xEy_row(Norm_A[i], Norm_B[i], &A[i * 9]);

	int iCount, iResult;
	iResult = Solve_Linear_Solution_Construction<_T>(A, 8, 9, NULL, x, &iCount);
	if (!iResult || iCount != 1)
		return 0;

	//分解
	_T U[3 * 3], S[3], Vt[3 * 3];
	SVD_Decompose(x, 3, S, U, Vt);
	//重新构造F
	memset(F, 0, 3 * 3 * sizeof(_T));
	F[0] = S[0], F[4] = S[1];
	Matrix_Multiply_3x3(U, F, F);
	Matrix_Multiply_3x3(F, Vt, F);

	//反归一化
	_T K9[3 * 3];
	K3_2_K9(T3_B, K9);
	Matrix_Transpose(K9, 3, 3, Temp);
	Matrix_Multiply_3x3(Temp, F, F);
	K3_2_K9(T3_A, K9);
	Matrix_Multiply_3x3(F, K9, F);

	//最后哦归一化，防止数据太小
	Normalize(F, 9, F);
	//Disp(F, 3, 3, "F");
	return 1;
}

//template void Gen_H_Coeff_row(double P1[2], double P2[2], double row[2 * 9]);
template<typename _T>void Gen_H_Coeff_row(_T P1[2], _T P2[2], _T row[2 * 9])
{//构造   uv2' x H * uv1 中的一行
	//此时，uv1, uv2 用的是齐次坐标，刚好等于归一化平面上的位姿看成空间点
	_T row_1[2 * 9] = { 0,0,0,-P1[0],-P1[1],-1,P1[0] * P2[1],P1[1] * P2[1],P2[1],
		P1[0],P1[1],1,0,0,0,-P1[0] * P2[0],-P1[1] * P2[0],-P2[0] };
	memcpy(row, row_1, 2 * 9 * sizeof(_T));
	return;
}

template<typename _T>_T Test_H(_T H[3 * 3], _T Point_A[2], _T Point_B[2])
{
	//验算 x2 x H * x1，理论上数值应该很小
	_T A[3] = { Point_A[0], Point_A[1], 1 },
		B[3] = { Point_B[0], Point_B[1], 1 },
		C[3];
	Matrix_Multiply_3x1(H, A, C);
	Cross_Product(B, C, C);

	//printf("x2'Hx1=%e\n", fError);

	Matrix_Multiply_3x1(H, A, C);
	Vector_Multiply(C, 3, 1 / C[2], C);
	_T fError = fGet_Distance(B, C, 3);
	//printf("x2 s= H*x1:%e\n", fError);

	return fError;
}
template<typename _T>_T Test_H(_T H[3 * 3], _T Point_A[][2], _T Point_B[][2], int iCount = 4)
{
	_T fError_Sum = 0;
	for (int i = 0; i < iCount; i++)
		fError_Sum += Test_H(H, Point_A[i], Point_B[i]);
	return fError_Sum;
}

//template void H_2_R_t(double H[3 * 3], double R[4][3 * 3], double t[4][3], double n[4][3]);
template<typename _T>void H_2_R_t(_T H[3 * 3], _T R[4][3 * 3], _T t[4][3], _T n[4][3])
{//理论上，一个H可以分解出4组R,t，还得进一步验最优解
	//首先，对H进行svd分解
	_T U[3 * 3], S[3], Vt[3 * 3];
	SVD_Decompose(H, 3, S, U, Vt);

	//分离出lambda
	_T lambda[3], lambda_sqr[3];
	Vector_Multiply(S, 3, 1 / S[1], lambda);
	Vector_Multiply(lambda, 3, lambda, lambda_sqr);

	_T n1, n3;
	_T x1, x3;

	//n1^2 = (λ1^2 -1)/(λ1^2 - λ3^2)
	n1 = (lambda_sqr[0] - 1) / (lambda_sqr[0] - lambda_sqr[2]);
	n3 = (1 - lambda_sqr[2]) / (lambda_sqr[0] - lambda_sqr[2]);
	x1 = sqrt(n1);
	x3 = sqrt(n3);

	//printf("x1:%e x3:%f\n", x1, x3);
	n[0][1] = n[1][1] = n[2][1] = n[3][1] = 0;
	n[0][0] = x1, n[0][2] = x3;
	n[1][0] = x1, n[1][2] = -x3;
	n[2][0] = -x1, n[2][2] = x3;
	n[3][0] = -x1, n[3][2] = -x3;
	//Disp((_T*)n, 4, 3, "n");

	_T lambda1_minus_lambda3 = lambda[0] - lambda[2];
	t[0][1] = t[1][1] = t[2][1] = t[3][1] = 0;
	t[0][0] = lambda1_minus_lambda3 * x1, t[0][2] = -lambda1_minus_lambda3 * x3;
	t[1][0] = t[0][0], t[1][2] = -t[0][2];
	t[2][0] = -t[0][0], t[2][2] = t[0][2];
	t[3][0] = -t[0][0], t[3][2] = -t[0][2];
	//Disp((_T*)t, 4, 3, "t");

	//构造R，0,3一样，2，4一样
	R[0][1] = R[1][1] = R[2][1] = R[3][1] =
		R[0][3] = R[1][3] = R[2][3] = R[3][3] =
		R[0][5] = R[1][5] = R[2][5] = R[3][5] =
		R[0][7] = R[1][7] = R[2][7] = R[3][7] = 0;
	R[0][4] = R[1][4] = R[2][4] = R[3][4] = 1;
	R[0][0] = R[1][0] = R[2][0] = R[3][0] = lambda[0] - lambda1_minus_lambda3 * x1 * x1;
	R[0][8] = R[1][8] = R[2][8] = R[3][8] = lambda[2] + lambda1_minus_lambda3 * x3 * x3;

	R[0][2] = -x1 * x3 * lambda1_minus_lambda3;
	R[0][6] = -R[0][2];

	R[1][2] = -R[0][2];
	R[1][6] = R[0][2];

	R[2][2] = -R[0][2];
	R[2][6] = R[0][2];

	R[3][2] = R[0][2];
	R[3][6] = -R[0][2];

	//for (int i = 0; i < 4; i++)
	//{
	//    Disp(R[i], 3, 3, "R");
	//    //Disp(n[i], 1, 3, "n");
	//    //Disp(t[i], 1, 3, "t");
	//}

	//回代R = U * Rsvd *  V'
	for (int i = 0; i < 4; i++)
	{
		Matrix_Multiply_3x3(U, R[i], R[i]);
		Matrix_Multiply_3x3(R[i], Vt, R[i]);

		Matrix_Multiply_3x1(U, t[i], t[i]);
		Matrix_Multiply(n[i], 1, 3, Vt, 3, n[i]);
		/*Disp(R[i], 3, 3, "R");
		Disp(t[i], 1, 3, "t");
		Disp(n[i], 1, 3, "n");*/
	}

	return;
}

template int Elect_R_t(double Point_A[5][2], double Point_B[5][2], int iPoint_Count, double R[][3 * 3], double t[][3], int iRt_Count, int* piBest_E_Index, double* pfError);
template<typename _T>int Elect_R_t(_T Point_A[5][2], _T Point_B[5][2], int iPoint_Count,
	_T R[][3 * 3], _T t[][3], int iRt_Count, int* piBest_E_Index,_T *pfError)
{
	int bFound = 0, iBest_E_Index = -1;
	_T fError_Sum, fBest_Error = MAX_FLOAT;
	for (int k = 0; k < iRt_Count; k++)
	{
		int bInvalid = 0;
		fError_Sum = 0;

		_T T[3 * 4], P[4];
		Gen_Pose_By_R_t(R[k], t[k], T, 0);
		/*if(k==3)
		Disp(T, 3, 4, "T");*/
		for (int j = 0; j < iPoint_Count; j++)
		{
			Triangulate_Cramer<_T>(Point_A[j], Point_B[j], T, P);
			Triangulate_Gauss<_T>(Point_A[j], Point_B[j], NULL, T, P);
			//if (k == 3 && j == 2)
				//Disp(P, 3, 1, "P");
			if (P[2] <= 0)
			{
				bInvalid = 1;
				break;
			}
			_T P2[3];
			TP(T, P, P2);
			fError_Sum += Test_Triangulate(P, T, Point_A[j], Point_B[j]);
			if (P2[2] <= 0)
			{
				bInvalid = 1;
				break;
			}
		}
		if (!bInvalid)
		{
			bFound++;
			if (fError_Sum < fBest_Error)
			{//成功了
				fBest_Error = fError_Sum;
				*piBest_E_Index = k;
			}
		}
	}

	if (pfError)
		*pfError = fBest_Error;
	if (fBest_Error == MAX_FLOAT)
		bFound = 0;
	return !!bFound;
}

template int Estimate_H_4Point(double Point_A[4][2], double Point_B[4][2], double T[3 * 4], double H[3 * 3], double* pfError);
template<typename _T>int Estimate_H_4Point(_T Point_A[4][2], _T Point_B[4][2], _T T[3 * 4], _T H[3 * 3],_T *pfError)
{//一般空间平面点求解H矩阵，理论上既可以算归一化平面，也可以算像素平面的点
	union {
		_T A[8 * 9];
		_T x[8 * 9];
	};
	for (int i = 0; i < 4; i++)
		Gen_H_Coeff_row(Point_A[i], Point_B[i], &A[i * 2 * 9]);

	//Disp(A, 8, 9, "A");
	int iCount, iResult;
	iResult = Solve_Linear_Solution_Construction<_T>(A, 8, 9, NULL, x, &iCount);
	if (!iResult || iCount != 1)
		return 0;   //病态数据，拒绝

	////H 矩阵归一化，非必要
	//for (int i = 0; i < 9;i++)
		//x[i] /= x[8];
	//printf("Error:%e\n", Test_H(x, Point_A, Point_B, 4));

	//注意，此处的t 也是一个scale，但是用法跟E矩阵的t不一样
	//此处的|t1|!= 1。 t1 = t/d，t 才是真实的位移
	//如果两个相机之间的位移为2，那么真实位移为 t = t1/|t1| * 2
	_T R[4][3 * 3], t[4][3], n[4][3];
	H_2_R_t(x, R, t, n);
	

	//接着还要三角化验算一把
	//注意，三角化的P点位置也只是尺度不确定未知，
	//要量出两相机之间的距离s 后，P = P * s/|t1| 才是
	//真实位置
	int iBest_Index;
	if (!(iResult = Elect_R_t(Point_A, Point_B, 4, R, t, 4, &iBest_Index,pfError)))
		return 0;

	if (T)
		Gen_Pose_By_R_t(R[iBest_Index], t[iBest_Index], T, 0);

	if (H)
	{   //H' 	= R + t1 * n'
		Matrix_Multiply(t[iBest_Index], 3, 1, n[iBest_Index], 3, H);
		Vector_Add(R[iBest_Index], H, 3 * 3, H);

		//H 矩阵归一化，非必要
		Vector_Multiply(H, 9, 1 / H[8], H);

		//printf("Error:%e\n", Test_H(H, Point_A, Point_B, 4));
		//Disp(H, 3, 3, "H");
	}
	return 1;
}

template double Sampson(double x1[2], double x2[2], double E[3 * 3]);
template<typename _T>_T Sampson(_T x1[2], _T x2[2], _T E[3 * 3])
{//算一个sampson 距离 (x2'*E*x1)^2/(Ex1[0]^2 +  Ex1[1]^2 + Ex2[0]^2 + Ex2[1]^2)
	_T Ex1[3], Ex2[3];
	//E * x1
	Ex1[0] = E[0] * x1[0] + E[1] * x1[1] + E[2];
	Ex1[1] = E[3] * x1[0] + E[4] * x1[1] + E[5];
	Ex1[2] = E[6] * x1[0] + E[7] * x1[1] + E[8];

	//E'*x2
	Ex2[0] = E[0] * x2[0] + E[3] * x2[1] + E[6];
	Ex2[1] = E[1] * x2[0] + E[4] * x2[1] + E[7];
	//Ex2[2] = E[2] * x2[0] + E[5] * x2[1] + E[8];

	_T x2_E_x1 = x2[0] * Ex1[0] + x2[1] * Ex1[1] + Ex1[2];
	return  x2_E_x1 * x2_E_x1 / (Ex1[0] * Ex1[0] +
		Ex1[1] * Ex1[1] + Ex2[0] * Ex2[0] + Ex2[1] * Ex2[1]);
}

template<typename _T>_T Sampson(_T x1[][2], _T x2[][2], int iCount, _T E[3 * 3])
{
	_T fError = 0;
	for (int i = 0; i < iCount; i++)
		fError += Sampson(x1[i], x2[i], E);
	return fError;
}
template void Sample_XY(double Point_A[][2], double Point_B[][2], int iPoint_Count, int iSample_Count);
template<typename _T>void Sample_XY(_T Point_A[][2], _T Point_B[][2], int iPoint_Count, int iSample_Count)
{//随机升起点对
	//出去以后，头iSample_Count 个点对必然是随机生成点对

	for (int i = 0; i < iSample_Count; i++)
	{//不必使用Colmap的算法
		int iRandom_Pos = iRandom() % (iPoint_Count - i) + i;
		//printf("%d\n", iRandom_Pos);
		//swap Pos i with iRandom_Pos
		SWAP(_T, Point_A[i][0], Point_A[iRandom_Pos][0]);
		SWAP(_T, Point_A[i][1], Point_A[iRandom_Pos][1]);
		SWAP(_T, Point_B[i][0], Point_B[iRandom_Pos][0]);
		SWAP(_T, Point_B[i][1], Point_B[iRandom_Pos][1]);
	}
	return;
}

template<typename _T>_T Compute_Squared_Sampson_Error(_T M[3 * 3], _T Point_A[][2], _T Point_B[][2], int iCount, _T eps,
	int* piInlier_Count, int bSwap_Inlier_Forward = 1)
{//计算所有点的
	int i, j = 0;
	_T fSum = 0;
	_T fError, fError_Sum = 0;

	//似乎没啥用
	//int iSeq_Not_Match_Count = 0;
	//const int Max_Seq_Nor_Match_Count = 50;
	for (i = j = 0; i < iCount /*&& iSeq_Not_Match_Count<= Max_Seq_Nor_Match_Count*/; i++)
	{
		if ((fError = Sampson(Point_A[i], Point_B[i], M)) <= eps)
		{//误差足够小，换到前面
			if (bSwap_Inlier_Forward)
			{
				SWAP(_T, Point_A[i][0], Point_A[j][0]);
				SWAP(_T, Point_A[i][1], Point_A[j][1]);
				SWAP(_T, Point_B[i][0], Point_B[j][0]);
				SWAP(_T, Point_B[i][1], Point_B[j][1]);
			}
			j++;
			fError_Sum += fError;
			//iSeq_Not_Match_Count = 0;
		}
		/*else
			iSeq_Not_Match_Count++;*/
	}

	*piInlier_Count = j;
	return fError_Sum;
}

int Ransac_Remain_Count(int iInlier, int iSample_Count, int iMin_Sample_Count)
{
	//内殿率
	float fInlier_Ration = (float)iInlier / iSample_Count;
	//n点全是内殿的概率
	float fP_Good = (float)pow(fInlier_Ration, iMin_Sample_Count);
	//n点不能构成E矩阵的概率，相当于Colmap 的denom
	float fP_Bad = 1 - fP_Good;
	//float fP_Fail_All = pow(fP_Bad, k);
	float fP = 0.99f;        //置信度

	if (fP_Bad <= 0)
		return 0;
	if (fP_Bad == 1.0)	//此时表示Inlier_Ratio=0，找不到任何匹配点
		return 0xFFFFFFF;	//返回个最大值

	int k = (int)ceil(3.f * log(1 - fP) / log(fP_Bad));
	return k;
}

template int Get_Inlier_Count(double Point_A[][2], double Point_B[][2], int iCount, double T[3 * 4], int bSwap_Inlier_Forward, int iMethod, double* pfError, double eps);
template<typename _T>int Get_Inlier_Count(_T Point_A[][2], _T Point_B[][2], int iCount, _T T[3 * 4], int bSwap_Inlier_Forward, int iMethod, _T* pfError, _T eps)
{//对一堆点对就T进行三角化，算出内殿数量
//三角化方法：0：DLT        SVD法
//            1: Cramer     克莱姆法
//            2: Gauss      迭代法

	_T fError, fError_Sum = 0;
	int iInlier, j;
	for (j = 0, iInlier = 0; j < iCount; j++)
	{
		/*if ( abs(Point_B[j][0]- -0.787080)<0.0001)
			printf("here");*/

		_T P[3], P2[3];
		if (iMethod == 0)
			Triangulate_DLT<_T>(Point_A[j], Point_B[j], NULL, T, P);
		else
		{
			Triangulate_Cramer<_T>(Point_A[j], Point_B[j], T, P);
			if (iMethod == 2)
				Triangulate_Gauss<_T>(Point_A[j], Point_B[j], NULL, T, P);
		}

		if (P[2] <= 0)
			continue;

		TP(T, P, P2);
		if (P2[2] <= 0)
			continue;

		fError = Test_Triangulate(P, T, Point_A[j], Point_B[j]);
		if (fError > eps)
			continue;
		//printf("Error:%f\n", fError);
		fError_Sum += fError;
		if (bSwap_Inlier_Forward)
		{
			//SWAP 到前面
			SWAP(_T, Point_A[j][0], Point_A[iInlier][0]);
			SWAP(_T, Point_A[j][1], Point_A[iInlier][1]);
			SWAP(_T, Point_B[j][0], Point_B[iInlier][0]);
			SWAP(_T, Point_B[j][1], Point_B[iInlier][1]);
		}
		iInlier++;
	}
	if (pfError)
		*pfError = fError_Sum;
	return iInlier;
}
template int Ransac_E(double Point_A[][2], double Point_B[][2], int iCount, double T[3 * 4], double E[3 * 3], int bUse_5_Point, int bPoint_In_Place);
template<typename _T>int Ransac_E(_T Point_A[][2], _T Point_B[][2], int iCount,
	_T T[3 * 4], _T E[3 * 3], int bUse_5_Point, int bPoint_In_Place)
{
	//为了不动原来点对，此处要赋值点对
	_T(*pDup_Point_A)[2], (*pDup_Point_B)[2];
	if (bPoint_In_Place)
	{
		pDup_Point_A = Point_A;
		pDup_Point_B = Point_B;
	}
	else
	{
		pDup_Point_A = (_T(*)[2])pMalloc(iCount * 2 * 2 * sizeof(_T));
		pDup_Point_B = pDup_Point_A + iCount;
		memcpy(pDup_Point_A, Point_A, iCount * 2 * sizeof(_T));
		memcpy(pDup_Point_B, Point_B, iCount * 2 * sizeof(_T));
	}

	const int iMax_Seq_Dup_Trail = 20, iHalf_Count = (iCount + 1) >> 1;
	const _T eps = 0.005208333333333333f * 0.005208333333333333f;

	int iIter, iRemain_Iter_Count = 0xFFFFFF, iSeq_Dup_Max_Inlier = 0;

	int iResult, iInlier, iMax_Inlier = 0, iSample_Count = bUse_5_Point ? 5 : 8;
	_T fError, fMin_Error = MAX_FLOAT;
	_T Best_E[3 * 3], T1[3 * 4], E1[3 * 3];

	//******************第一步，粗估计**********************************************************/
	for (iIter = 0; iIter < iRemain_Iter_Count && 
		iMax_Inlier< iHalf_Count &&			//良性数据能加快跳出
		iSeq_Dup_Max_Inlier < iMax_Seq_Dup_Trail; iIter++)
	{
		//随机选出一组点对
		Sample_XY(pDup_Point_A, pDup_Point_B, iCount, iSample_Count);

		if (bUse_5_Point)
		{
			if (!(iResult = Estimate_E_T_5Point(pDup_Point_A, pDup_Point_B, T1, E1)))
				continue;
		}else
		{
			if (!(iResult = Estimate_E_8Point<_T>(pDup_Point_A, pDup_Point_B, T1, E1)))
				continue;
		}

		fError = Compute_Squared_Sampson_Error(E1, pDup_Point_A, pDup_Point_B, iCount, eps, &iInlier, 0);
		if (iInlier < iSample_Count)
			continue;

		int iPre_Max_Inlier = iMax_Inlier;
		if ((iInlier > iMax_Inlier) ||
			(iInlier == iMax_Inlier && fError < fMin_Error))
		{
			memcpy(Best_E, E1, 3 * 3 * sizeof(_T));
			fMin_Error = fError;
			iMax_Inlier = iInlier;
			iSeq_Dup_Max_Inlier = 0;
			//fError = Compute_Squared_Sampson_Error(Best_E, pDup_Point_A, pDup_Point_B, iCount, eps, &iInlier, 0);
		}
		else if (iPre_Max_Inlier == iMax_Inlier)    //若然连续多次冲不破最大内点数，则跳出
			iSeq_Dup_Max_Inlier++;
				
		//printf("Inlier:%d\n", iPre_Max_Inlier);
		iRemain_Iter_Count = Ransac_Remain_Count(iMax_Inlier, iCount, iSample_Count);
	}
	//******************第一步，粗估计**********************************************************/

	//升起全部Inlier
	fError = Compute_Squared_Sampson_Error(Best_E, pDup_Point_A, pDup_Point_B, iCount, eps, &iInlier);

	//******************************第二步，精估计，步步三角化***********************************/
	iMax_Inlier = 0;
	fMin_Error = MAX_FLOAT;
	for (int i = 0; i < 5; i++)
	{
		if (!(iResult = Estimate_E_T_nPoint<_T>(pDup_Point_A, pDup_Point_B, iInlier, &iInlier, T1, E1)))
			break;
		iInlier = Get_Inlier_Count(pDup_Point_A, pDup_Point_B, iCount, T1, 1, 2, &fError);
		printf("Inlier:%d average:%f\n", iInlier, fError / iInlier);
		if (iMax_Inlier < iInlier ||
			(iMax_Inlier == iInlier && fMin_Error < fError))
		{
			iMax_Inlier = iInlier;
			fMin_Error = fError;
			if (T)
				memcpy(T, T1, 3 * 4 * sizeof(_T));
			if (E)
				memcpy(E, E1, 3 * 3 * sizeof(_T));
		}else
			break;
	}
	//******************************第二步，精估计，步步三角化***********************************/

	if (!bPoint_In_Place)
		Free(pDup_Point_A);
	return iMax_Inlier >= iSample_Count ? 1 : 0;
}

template<typename _T>_T Computer_H_Error(_T H[3 * 3], _T Point_A[2], _T Point_B[2])
{//返回经过H变换后的欧几里得距离（误差）eps 为阀值
	_T d0 = H[0 * 3 + 0] * Point_A[0] + H[0 * 3 + 1] * Point_A[1] + H[0 * 3 + 2],
		d1 = H[1 * 3 + 0] * Point_A[0] + H[1 * 3 + 1] * Point_A[1] + H[1 * 3 + 2],
		d2 = H[2 * 3 + 0] * Point_A[0] + H[2 * 3 + 1] * Point_A[1] + H[2 * 3 + 2];

	d0 = d0 / d2 - Point_B[0];
	d1 = d1 / d2 - Point_B[1];
	return d0 * d0 + d1 * d1;
}

template<typename _T>_T Computer_H_Error(_T H[3 * 3], _T Point_A[][2], _T Point_B[][2], int iCount, _T eps,
	int* piInlier_Count = NULL, int bSwap_Inlier_Forward = 1)
{//测试一群点对经过H 变换后的中体误差
	_T fError, fError_Sum = 0;
	int j = 0;
	for (int i = 0; i < iCount; i++)
	{
		/*if ( abs(Point_B[i][0]- 0.388640f)<0.001 ||
			abs(Point_B[i][0] - 0.812960f) < 0.001 ||
			abs(Point_B[i][0] - -0.042000) < 0.001 ||
			abs(Point_B[i][0] - 0.380480f) < 0.001 )
		{
			printf("%f %f %f %f\n", Point_A[i][0], Point_A[i][1], Point_B[i][0], Point_B[i][1]);
		}*/

		fError = Computer_H_Error(H, Point_A[i], Point_B[i]);
		//printf("i:%d Error:%f\n", i, fError);
		if (fError < eps)
		{//及格，将点移前
			if (bSwap_Inlier_Forward)
			{
				SWAP(_T, Point_A[i][0], Point_A[j][0]);
				SWAP(_T, Point_A[i][1], Point_A[j][1]);
				SWAP(_T, Point_B[i][0], Point_B[j][0]);
				SWAP(_T, Point_B[i][1], Point_B[j][1]);
			}
			j++;
			fError_Sum += fError;
		}
	}

	if (piInlier_Count)
		*piInlier_Count = j;

	return fError_Sum;
}

template<typename _T>int Estimate_H_nPoint(_T Point_A[5][2], _T Point_B[5][2], int iCount, int* piInlier,
	_T T[3 * 4] = NULL, _T H[3 * 3] = NULL, _T eps = 0.0002f)
{
	int bRet = 0, iInlier, iMax_Inlier = 0, iBest_Index = -1;
	_T fMin_Error = MAX_FLOAT, fError_Sum, * pA = (_T*)pMalloc(iCount * 2 * 9 * sizeof(_T));
	if (!pA || iCount < 4)
		goto END;

	for (int i = 0; i < iCount; i++)
	{
		Gen_H_Coeff_row(Point_A[i], Point_B[i], &pA[i * (2 * 9)]);
		//printf("%f %f %f %f\n", Point_A[i][0], Point_A[i][1], Point_B[i][0], Point_B[i][1]);
	}

	int iResult;
	union {
		_T x[9];
		_T R1[3 * 3];
	};
	Solve_Linear_Contradictory<_T>(pA, iCount * 2, 9, NULL, x, &iResult);
	if (!iResult)
		return 0;

	//接着需要分解R,t
	_T(*T2)[3 * 4], (*n)[3];    //此处应该借空间
	T2 = (_T(*)[3 * 4])pA;
	n = (_T(*)[3]) (T2 + 4);
	{
		_T R[4][3 * 3], t[4][3];
		H_2_R_t(x, R, t, n);
		Gen_Pose_By_R_t<_T>(R[0], t[0], T2[0], 0);
		Gen_Pose_By_R_t(R[1], t[1], T2[1], 0);
		Gen_Pose_By_R_t(R[2], t[2], T2[2], 0);
		Gen_Pose_By_R_t(R[3], t[3], T2[3], 0);
	}

	//接着还要三角化验算一把
	//注意，三角化的P点位置也只是尺度不确定未知，
	//要量出两相机之间的距离s 后，P = P * s/|t1| 才是
	//真实位置
	for (int k = 0; k < 4; k++)
	{
		fError_Sum = 0; iInlier = 0;
		_T* pT = T2[k], P[4];
		for (int j = 0; j < iCount; j++)
		{
			Triangulate_Cramer<_T>(Point_A[j], Point_B[j], pT, P);
			Triangulate_Gauss<_T>(Point_A[j], Point_B[j], NULL, pT, P);;
			if (P[2] <= 0)
				continue;

			_T P2[3];
			TP(pT, P, P2);
			if (P2[2] <= 0)
				continue;

			_T  fError = Test_Triangulate(P, pT, Point_A[j], Point_B[j]);
			if (fError > eps)
				continue;

			fError_Sum += fError;
			iInlier++;
		}

		if (iInlier > iMax_Inlier ||
			(iInlier == iMax_Inlier && fError_Sum < fMin_Error))
		{
			iMax_Inlier = iInlier;
			iBest_Index = k;
			fMin_Error = fError_Sum;
		}
	}

	if (iBest_Index == -1)
		goto END;
	if (piInlier)
		*piInlier = iMax_Inlier;
	if (T)
		memcpy(T, T2[iBest_Index], 3 * 4 * sizeof(_T));

	if (H)
	{//需要生成一个新的H 矩阵
		_T* t1 = (_T*)(n + 4),
			* n1 = t1 + 3;
		memcpy(n1, n[iBest_Index], 3 * sizeof(_T));
		Get_R_t(T2[iBest_Index], R1, t1);

		//H' 	= R + t1 * n'
		Matrix_Multiply(t1, 3, 1, n1, 3, H);
		Vector_Add(R1, H, 3 * 3, H);
	}

	bRet = 1;
END:
	Free(pA);
	return bRet;
}
template int Ransac_H(double Point_A[][2], double Point_B[][2], int iCount, double T[3 * 4], double H[3 * 3], int bPoint_In_Place, double eps);
template<typename _T>int Ransac_H(_T Point_A[][2], _T Point_B[][2], int iCount,	_T T[3 * 4], _T H[3 * 3], int bPoint_In_Place, _T eps)
{
	//为了不动原来点对，此处要赋值点对
	_T(*pDup_Point_A)[2], (*pDup_Point_B)[2];
	if (bPoint_In_Place)
	{
		pDup_Point_A = Point_A;
		pDup_Point_B = Point_B;
	}
	else
	{
		pDup_Point_A = (_T(*)[2])pMalloc(iCount * 2 * 2 * sizeof(_T));
		pDup_Point_B = pDup_Point_A + iCount;
		memcpy(pDup_Point_A, Point_A, iCount * 2 * sizeof(_T));
		memcpy(pDup_Point_B, Point_B, iCount * 2 * sizeof(_T));
	}

	const int iMax_Seq_Dup_Trail = 20, iHalf_Count = (iCount + 1) >> 1;
	//const _T eps = 0.0004;   // 0.005208333333333333f * 0.005208333333333333f;

	int iIter, iRemain_Iter_Count = 0xFFFFFF, iSeq_Dup_Max_Inlier = 0;

	int iResult, iInlier, iMax_Inlier = 0, iSample_Count = 4;
	_T fError, fMin_Error = MAX_FLOAT;
	_T Best_H[3 * 3], T1[3 * 4], H1[3 * 3];

	//******************第一步，粗估计**********************************************************/
	for (iIter = 0; iIter < iRemain_Iter_Count
		&& iMax_Inlier < iHalf_Count   			//良性数据能加快跳出
		&& iSeq_Dup_Max_Inlier < iMax_Seq_Dup_Trail; iIter++)
	{
		//随机选出一组点对
		Sample_XY(pDup_Point_A, pDup_Point_B, iCount, iSample_Count);
		_T T[3 * 4];

		if (!(iResult = Estimate_H_4Point(pDup_Point_A, pDup_Point_B, T, H1)))
			continue;

		fError = Computer_H_Error(H1, Point_A, Point_B, iCount, eps, &iInlier);

		if (iInlier < iSample_Count)
			continue;

		int iPre_Max_Inlier = iMax_Inlier;
		if ((iInlier > iMax_Inlier) ||
			(iInlier == iMax_Inlier && fError < fMin_Error))
		{
			memcpy(Best_H, H1, 3 * 3 * sizeof(_T));
			memcpy(T1, T, 3 * 4 * sizeof(_T));
			fMin_Error = fError;
			iMax_Inlier = iInlier;
			iSeq_Dup_Max_Inlier = 0;
		}
		else if (iPre_Max_Inlier == iMax_Inlier)    //若然连续多次冲不破最大内点数，则跳出
			iSeq_Dup_Max_Inlier++;

		//printf("Inlier:%d\n", iPre_Max_Inlier);
		iRemain_Iter_Count = Ransac_Remain_Count(iMax_Inlier, iCount, iSample_Count);
	}
	//******************第一步，粗估计**********************************************************/

	//升起全部Inlier
	fError = Computer_H_Error(Best_H, pDup_Point_A, pDup_Point_B, iCount, eps, &iInlier);

	//用全体三角化方法
	//iInlier = Get_Inlier_Count(pDup_Point_A, pDup_Point_B, iCount, T1, 1, 2, &fError, eps);
	//printf("Corse Estimate:%e\n", fError/iInlier);


	//******************************第二步，精估计，步步三角化***********************************/
	iMax_Inlier = -1;
	fMin_Error = MAX_FLOAT;
	for (int i = 0; i < 5; i++)
	{
		//用及格的点你和H,位姿
		if (!(iResult = Estimate_H_nPoint(Point_A, Point_B, iInlier, &iInlier, T1, H1, eps)))
			break;

		iInlier = Get_Inlier_Count(pDup_Point_A, pDup_Point_B, iCount, T1, 1, 2, &fError, eps);
		printf("Inlier:%d average:%e\n", iInlier, fError / iInlier);
		if (iMax_Inlier < iInlier ||
			(iMax_Inlier == iInlier && fError < fMin_Error))
		{
			iMax_Inlier = iInlier;
			fMin_Error = fError;
			if (T)
				memcpy(T, T1, 3 * 4 * sizeof(_T));
			if (H)
				memcpy(H, H1, 3 * 3 * sizeof(_T));
		}
		else
			break;
	}

	if (!bPoint_In_Place)
		Free(pDup_Point_A);
	return iMax_Inlier >= iSample_Count ? 1 : 0;
}
//*********************************求E矩阵**************************/

//*****************************三种三角化*********************************/
//这些小块展开算法应该还有利用价值
#define AtA_3x2(A,AtA)  \
{                   \
    AtA[0] = A[0] * A[0] + A[2] * A[2] + A[4] * A[4];   \
    AtA[1] = A[0] * A[1] + A[2] * A[3] + A[4] * A[5];   \
    AtA[2] = A[1] * A[0] + A[3] * A[2] + A[5] * A[4];   \
    AtA[3] = A[1] * A[1] + A[3] * A[3] + A[5] * A[5];   \
}\

#define Atb_2x3_2x1(A,b,Atb)  \
{\
    Atb[0] = A[0] * b[0] + A[2] * b[1] + A[4] * b[2];   \
    Atb[1] = A[1] * b[0] + A[3] * b[1] + A[5] * b[2];   \
}\

template int Triangulate_Cramer(double Pn1[3], double Pn2[3], double T[3 * 4], double P[3]);
template<typename _T>int Triangulate_Cramer(_T Pn1[3], _T Pn2[3], _T T[3 * 4], _T P[3])
{//就是一般方程法，除了快没别的好，通常用作初值
	//R * Pn1	Pn2	*	s1	=	-t      就是解这个方程
	//                  s2
	_T x[2];
	{//第一步，算出 s1,s2
		_T A[3 * 2], b[3];
		//R* Pn1
		A[0] = T[0] * Pn1[0] + T[1] * Pn1[1] + T[2];;
		A[2] = T[4] * Pn1[0] + T[5] * Pn1[1] + T[6];
		A[4] = T[8] * Pn1[0] + T[9] * Pn1[1] + T[10];

		//-Pn2
		A[1] = -Pn2[0], A[3] = -Pn2[1], A[5] = -1;

		//b = -t
		b[0] = -T[3], b[1] = -T[7], b[2] = -T[11];

		//此处用快速方法
		_T AtA[2 * 2], Atb[2];
		//A' =  0	2	4       A = 0   1           
		//      1	3	5           2   3
		//                          4   5
		AtA_3x2(A, AtA);
		Atb_2x3_2x1(A, b, Atb);

		if (!Solve_Linear_2x2(AtA, Atb, x))
			return 0;
	}

	_T P1[3], P2[4];
	//P1 = s1 * Pn1
	P1[0] = Pn1[0] * x[0], P1[1] = Pn1[1] * x[0], P1[2] = x[0];

	//P2 = s2 * Pn2     注意，P2 还是相机坐标
	P2[0] = Pn2[0] * x[1], P2[1] = Pn2[1] * x[1], P2[2] = x[1], P2[3] = 1;

	//P2 =>世界坐标
	_T T_Inv[4 * 4];
	Get_Inv_T(T, T_Inv);
	TP(T_Inv, P2, P2);

	P[0] = (P1[0] + P2[0]) * 0.5f;
	P[1] = (P1[1] + P2[1]) * 0.5f;
	P[2] = (P1[2] + P2[2]) * 0.5f;
	return 1;
}
template int Triangulate_DLT(double Pn1[3], double Pn2[3], double T1[3 * 4], double T2[3 * 4], double P[3]);
template<typename _T>int Triangulate_DLT(_T Pn1[3], _T Pn2[3], _T T1[3 * 4], _T T2[3 * 4], _T P[3])
{//给定两个归一化平面点，两个位姿，需按照空间点P
	_T A[4 * 4] = {};    //解一个矛盾方程组
	if (T1)
	{
		//相机1
		//xn*P3 - P1 *	X	= 0
		Vector_Multiply<_T>(&T1[2 * 4], 4, Pn1[0], &A[0 * 4]);
		Vector_Minus(&A[0 * 4], &T1[0 * 4], 4, &A[0 * 4]);

		//yn*P3 - P2
		Vector_Multiply(&T1[2 * 4], 4, Pn1[1], &A[1 * 4]);
		Vector_Minus(&A[1 * 4], &T1[1 * 4], 4, &A[1 * 4]);
	}
	else
	{//第1，2行可以快速赋值
		A[0] = A[5] = -1;
		A[2] = Pn1[0], A[6] = Pn1[1];
	}

	//相机2
	//xn*P3 - P1 *	X	= 0
	//Vector_Multiply<_T>(&T2[2 * 4], 4, Pn2[0], &A[2 * 4]);
	//Vector_Minus(&A[2 * 4], &T2[0 * 4], 4, &A[2 * 4]);
	//直接展开
	A[2 * 4 + 0] = T2[2 * 4 + 0] * Pn2[0] - T2[0 * 4 + 0];
	A[2 * 4 + 1] = T2[2 * 4 + 1] * Pn2[0] - T2[0 * 4 + 1];
	A[2 * 4 + 2] = T2[2 * 4 + 2] * Pn2[0] - T2[0 * 4 + 2];
	A[2 * 4 + 3] = T2[2 * 4 + 3] * Pn2[0] - T2[0 * 4 + 3];

	//yn*P3 - P2
	//Vector_Multiply(&T2[2 * 4], 4, Pn2[1], &A[3 * 4]);
	//Vector_Minus(&A[3 * 4], &T2[1 * 4], 4, &A[3 * 4]);
	//Disp(A, 4, 4, "A");
	A[3 * 4 + 0] = T2[2 * 4 + 0] * Pn2[1] - T2[1 * 4 + 0];
	A[3 * 4 + 1] = T2[2 * 4 + 1] * Pn2[1] - T2[1 * 4 + 1];
	A[3 * 4 + 2] = T2[2 * 4 + 2] * Pn2[1] - T2[1 * 4 + 2];
	A[3 * 4 + 3] = T2[2 * 4 + 3] * Pn2[1] - T2[1 * 4 + 3];
	//Disp(A, 4, 4, "A");

	_T X[4];
	int iResult;
	//注意，此处无计可施，方程太少，不能用反幂法
	//所以这是个龟慢的算法
	Solve_Homo_Linear_SVD(A, 4, 4, X, &iResult);
	Vector_Multiply(X, 3, 1.f / X[3], P);
	//Disp(P, 1, 4, "P");
	return 0;
}


//A' =	0	3	6	9		A = 0	1	2		
//		1	4	7	10			3	4	5
//		2	5	8	11			6	7	8
//								9	10	11

#define AtA_4x3(J,H)  \
{                   \
    H[0] = J[0] * J[0] + J[3] * J[3] + J[6] * J[6] + J[9] * J[9];	\
	H[1] = H[3] = J[0] * J[1] + J[3] * J[4] + J[6] * J[7] + J[9] * J[10];	\
	H[2] = H[6] = J[0] * J[2] + J[3] * J[5] + J[6] * J[8] + J[9] * J[11];	\
\
	H[4] = J[1] * J[1] + J[4] * J[4] + J[7] * J[7] + J[10] * J[10];	\
	H[5] = H[7] = J[1] * J[2] + J[4] * J[5] + J[7] * J[8] + J[10] * J[11];	\
\
	H[8] = J[2] * J[2] + J[5] * J[5] + J[8] * J[8] + J[11] * J[11];	\
}\

template<typename _T>static void Trianglate_Get_H_JtE(_T P1[3], _T Pn1[2], _T R1[3 * 3],
	_T P2[3], _T Pn2[2], _T R2[3 * 3], _T H[], _T JtE[], _T E[3])
{
	_T J[4 * 3];
	Get_dPn_dP(P1, R1, J);
	Get_dPn_dP(P2, R2, &J[2 * 3]);
	//Transpose_Multiply(J, 4, 3, H, 0);
	AtA_4x3(J, H);

	E[0] = Pn1[0] - P1[0] / P1[2],
		E[1] = Pn1[1] - P1[1] / P1[2],
		E[2] = Pn2[0] - P2[0] / P2[2],
		E[3] = Pn2[1] - P2[1] / P2[2];

	//At_x_B(J, 4, 3, E, 1, JtE);
	JtE[0] = J[0] * E[0] + J[3] * E[1] + J[6] * E[2] + J[9] * E[3];
	JtE[1] = J[1] * E[0] + J[4] * E[1] + J[7] * E[2] + J[10] * E[3];
	JtE[2] = J[2] * E[0] + J[5] * E[1] + J[8] * E[2] + J[11] * E[3];
	return;
}
template int Triangulate_Gauss(double Pn1[3], double Pn2[3], double T1[4 * 4], double T2[4 * 4], double P[3]);
template<typename _T>int Triangulate_Gauss(_T Pn1[3], _T Pn2[3], _T T1[4 * 4], _T T2[4 * 4], _T P[3])
{//迭代法，由于找不到数据证明它有用，暂时不优化了。留待日后再搞
	_T fError, fPre_Error = MAX_FLOAT, fError_Delta, eps = 1e-10;
	const int iMax_Iter = 30;
	int bRet = 0, iIter, iResult;

	//设置求解参数初值
	_T P1[3], P2[3], P_Dup[3] = { P[0],P[1],P[2] },
		R1[3 * 3], R2[3 * 3];

	if (!T1)//如果将相机1视为世界坐标，则直接赋值
		Gen_I_Matrix(R1, 3, 3);
	else
		Get_R_t<_T>(T1, R1, NULL);
	Get_R_t<_T>(T2, R2, NULL);

	//P1[0] = P[0], P1[1] = P[1], P1[2] = P[2];
	for (iIter = 0; iIter < iMax_Iter; iIter++)
	{//迭代体
		if (T1)
			TP(T1, P_Dup, P1);
		else
			P1[0] = P_Dup[0], P1[1] = P_Dup[1], P1[2] = P_Dup[2];

		TP(T2, P_Dup, P2);
		//注意：以下将Pd 作为Pn 的初值
		_T H[3 * 3], JtE[3], E[4];
		Trianglate_Get_H_JtE(P1, Pn1, R1, P2, Pn2, R2, H, JtE, E);

		//误差检测
		fError = fGet_Sqr_Sum(E, 4);
		fError_Delta = fPre_Error - fError;
		if (fError < eps)
			break;
		if (fError_Delta < 0)
			break;   //发散了
		fPre_Error = fError;
		P[0] = P_Dup[0], P[1] = P_Dup[1], P[2] = P_Dup[2];

		_T Delta_x[3];
		_T A[4 + 3 + 2] = { H[0], H[1], H[2], JtE[0],
								H[4],H[5],JtE[1],
									H[8],JtE[2] };
		if (!(iResult = Solve_Linear_Gause_AAt<_T>(A, 3, NULL, Delta_x)))
			break;

		//更新P
		Vector_Add(P_Dup, Delta_x, 3, P_Dup);
	}

	return 1;
}

template void Test_E(double E[3 * 3], double P1[5][2], double P2[5][2], int iCount);
template<typename _T> void Test_E(_T E[3 * 3], _T P1[5][2], _T P2[5][2], int iCount)
{
	const _T eps = 1e-10;
	_T Temp[3], fTotal = 0;
	for (int i = 0; i < iCount; i++)
	{
		_T P_1[3] = { P1[i][0],P1[i][1],1 },
			P_2[3] = { P2[i][0], P2[i][1], 1 }, fValue;

		/*Matrix_Multiply(P_1, 1, 3, E, 3, Temp);
		fValue = fDot(Temp, P_2, 2);*/

		//改为 x2' * E * x1 = 0
		Matrix_Multiply(P_2, 1, 3, E, 3, Temp);
		fValue = fDot(Temp, P_1, 2);

		fTotal += fValue * fValue;
	}
	printf("第一步，验算x1'Ex2 的总体误差：Error:%e\n", fTotal);
	if (fTotal > eps)
		printf("过不了总体误差\n");

	//第二步，行列式约束：
	_T fValue = fGet_Determinant(E, 3);
	printf("行列式：%e\n", fValue);
	if (fValue > eps)
		printf("过不了行列式约束\n");

	//第三步，迹约束 2EE'E - trace(EE')E = 0
	_T EEtE_x_2[3 * 3];
	Transpose_Multiply(E, 3, 3, EEtE_x_2);
	_T fTrace = fGet_Tr(EEtE_x_2, 3);

	Matrix_Multiply_3x3(EEtE_x_2, E, EEtE_x_2);
	Vector_Multiply<_T>(EEtE_x_2, 3 * 3, 2, EEtE_x_2);

	union {
		_T Trace_E[3 * 3];
		_T Result[3 * 3];
	};
	Vector_Multiply(E, 3 * 3, fTrace, Trace_E);
	Vector_Minus(EEtE_x_2, Trace_E, 3 * 3, Result);
	Disp(Result, 3, 3, "2EE'E - trace(EE')E");
	if (fGet_Mod(Result, 9) > eps)
		printf("过不了迹约束:%e\n", fGet_Mod(Result, 9));
	return;
}

template double Test_Triangulate(double P[3], double T[4 * 4], double Point_A[2], double Point_B[2]);
template<typename _T>_T Test_Triangulate(_T P[3], _T T[4 * 4], _T Point_A[2], _T Point_B[2])
{//在归一化平面上检验三角化结果
	_T Pn1[2] = { P[0] / P[2], P[1] / P[2] };
	_T P2[4];
	TP(T, P, P2);
	_T Pn2[2] = { P2[0] / P2[2],P2[1] / P2[2] };
	_T fError = fGet_Distance(Pn1, Point_A, 2);
	fError += fGet_Distance(Pn2, Point_B, 2);
	/*printf("xn1:%f yn1:%f ref_x:%f ref_y:%f\nxn2:%f yn2:%f ref_x:%f ref_y:%f Error:%e\n",
		Pn1[0], Pn1[1], Ref_1[0], Ref_1[1],
		Pn2[0], Pn2[1], Ref_2[0], Ref_2[1],
		fError);*/
	return fError;
}

template double Test_Triangulate(double T[4 * 4], double Point_A[][2], double Point_B[][2], int iCount);
template<typename _T>_T Test_Triangulate(_T T[4 * 4], _T Point_A[][2], _T Point_B[][2], int iCount)
{//算一堆点对的三角化误差
	_T fTotal = 0, fError;
	_T R[3 * 3], t[3];
	Get_R_t(T, R, t);

	for (int i = 0; i < iCount; i++)
	{
		_T P[3];
		Triangulate_Cramer(Point_A[i], Point_B[i], T, P);
		Triangulate_Gauss<_T>(Point_A[i], Point_B[i], NULL, T, P);
		//Triangulate_DLT<_T>(Point_A[i], Point_B[i], NULL, T, P);
		fError = Test_Triangulate<_T>(P, T, Point_A[i], Point_B[i]);
		fTotal += fError;
	}
	//printf("Error:%e\n", fTotal);
	return fTotal;
}
//*****************************三种三角化*********************************/
