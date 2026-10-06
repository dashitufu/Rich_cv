//这组试验围绕着选手逆转展开
#include <cmath>
#include "Slam.h"
#include "Image.h"
#include "Chess_Board_Detect.h"

void Rotation_Test_1()
{//先看看能转多少度，哪只手系
	typedef double _T;
	Point_Cloud<_T> oPC;
	Init_Point_Cloud(&oPC, 10000);
	{
		const int iStep_Count = 100;
		_T Rotation_Vector[4] = { 1,1,1 }, P0[3] = { 0.5, 0, 0 };

		Normalize(Rotation_Vector, 3, Rotation_Vector);
		Draw_Line<_T>(&oPC, -1, -1,-1, 1, 1,1);
		

		for (int i = 0; i < iStep_Count; i++)
		{
			_T fAngle = (2 * PI) * i / iStep_Count;
			Rotation_Vector[3] = fAngle;
			_T R[3 * 3];
			Rotation_Vector_4_2_Matrix(Rotation_Vector, R);
			_T P1[3];
			Matrix_Multiply_3x1(R, P0, P1);
			Draw_Point(&oPC, P1[0], P1[1], P1[2], (int)(_T)i*255/iStep_Count,0,0 );
			Disp(Rotation_Vector, 1, 3, "R V");
		}
		bSave_PLY("c:\\tmp\\1.ply", oPC);
	}
	Free_Point_Cloud(&oPC);
}
void T_Inv_Test()
{//快速求T的逆
	typedef double _T;
	_T V[3] = { 1,2,3 }, t[3] = { 2,3,4 }, T[4 * 4], T_Inv[4 * 4];
	int iResult = 0;
	Gen_Pose_By_V3_t(V, t,T);
	Disp(T, 4, 4, "T");

	Get_Inv_Matrix(T, T_Inv, 4, &iResult);
	Disp(T_Inv, 4, 4, "T_Inv");
		
	_T R[3 * 3], t1[3];
	Vector_Multiply<_T>(t, 3, -1, t1);
	Rotation_Vector_3_2_Matrix(V, R);
	Matrix_Transpose(R, 3, 3, R);
	Matrix_Multiply_3x1(R, t1, t1);
	Gen_Pose_By_R_t(R, t1, T_Inv);
	Disp(T_Inv, 4, 4, "T_Inv");


	Get_Inv_T(T, T_Inv);
	Disp(T_Inv, 4, 4, "T_Inv");
	return;
}
void Draw_Pose_Test()
{
	typedef double _T;
	Point_Cloud<_T> oPC;
	Init_Point_Cloud(&oPC, 10000);
	_T V[4] = { 1,0,0,PI / 4 }, t[3] = { 2,3,4}, T[4 * 4];
	Gen_Pose_By_V3_t(V, t, T);
	Draw_Camera(&oPC, T);

	bSave_PLY("c:\\tmp\\1.ply", oPC);

	Free_Point_Cloud(&oPC);
}

template<typename _T>void Get_Jl_Inv_4(_T V[4], _T J_Inv[3 * 3])
{
	_T theta_div_2 = V[3]/2,
		cot_div_2 = theta_div_2 /tan(theta_div_2);
	union {
		_T aat[3 * 3];
		_T a_hat[3 * 3];
	};
		
	Transpose_Multiply(V, 3, 1, aat);
	_T A[3 * 3];
	Gen_I_Matrix(A, 3, 3);
	Vector_Multiply(A, 3 * 3, cot_div_2, A);

	Vector_Multiply(aat, 3 * 3, 1 - cot_div_2, aat);
	Vector_Add(A, aat, 3 * 3, A);

	Hat(V, a_hat);
	Vector_Multiply<_T>(a_hat,3*3,theta_div_2, a_hat);

	Vector_Minus(A, a_hat, 3 * 3, J_Inv);
	return;
}
template<typename _T>void Get_Jl_Inv_3(_T V[3], _T J_Inv[3 * 3])
{
	_T V1[4];
	Rotation_Vector_3_2_4(V, V1);
	Get_Jl_Inv_4(V1, J_Inv);
}
void Exp_Test_1()
{//尝试一下求Exp(A)
	typedef double _T;
	_T V[3] = { 0,0,PI*0.8 }, R[3 * 3];
	
	{
		Hat(V, R);
		Exp_M(R, 3, R, 1e-5);
		Disp(R, 3, 3, "R");

		Exp_V3(V, R);
		Disp(R, 3, 3, "R");
	}

	{
		Rotation_Vector_3_2_Matrix(V, R);
		Disp(R, 3, 3, "R");
	}

	{
		_T V1[3];
		Rotation_Matrix_2_Vector_3<_T>(R, V1);
		Disp(V1, 1, 3, "v");
	}
	{
		_T Ksi[6] = { 1,0,0,0,0,PI / 3 }, T[44*4];
		se3_2_SE3(Ksi, T);
		Disp(T, 4, 4, "T");
	}

	{//主要看J*rho 的轨迹。随R 前进的联动方式已经确立
		Point_Cloud<_T> oPC;
		Init_Point_Cloud(&oPC, 10000);
		_T Ksi[6] = { 0.8,0.8,0.8,0,0,PI / 3 }, T[4 * 4];
		for (int i = 0; i < 100; i++)
		{
			_T Ksi_1[6] = { 0.8*i/100,0.8 * i / 100,0.8 * i / 100,0,0,(PI*2)*i/100 };
			se3_2_SE3(Ksi_1, T);
			Draw_Point(&oPC, T[3], T[7], T[11]);
		}
		for (int i = 0; i < 100; i++)
			Draw_Point<_T>(&oPC, 0, 0, (_T)i / 100);
		//bSave_PLY("c:\\tmp\\1.ply", oPC);
		Free_Point_Cloud(&oPC);
	}

	{//BCH 上的雅可比，此处与矩阵求逆已经对齐
		_T V[3] = { 0.4,0.5,0.6 }, Jl[3 * 3];
		Get_Jl_3(V, Jl);
		Disp(Jl, 3, 3, "Jl");

		_T Jl_Inv[3 * 3];
		Get_Jl_Inv_3(V, Jl_Inv);
		Disp(Jl_Inv, 3, 3, "Jl_Inv");

		Get_Inv_Matrix(Jl, Jl_Inv, 3);
		Disp(Jl_Inv, 3, 3, "Jl_Inv");
	}

	{//delta_phi 进行左扰动，看 delta_phi + phi 是什么？
	//按照李代数的定义，本来在R^3 上是没有定义 delta_phi + phi 的
		_T phi[3] = { 1,2,3 },
			delta_phi[4] =  { 4, 3, 1, 0.05 },
			phi_1[3],  R[3 * 3], Rl[3 * 3], Temp[3 * 3];
		_T Jl[3 * 3];

		//首先讲delta 化为3维向量
		Normalize(delta_phi, 3, delta_phi);
		Rotation_Vector_4_2_3<_T>(delta_phi, delta_phi);
				
		//蜀商的公式，exp(delta_phi + phi) = exp[(J(phi)*delta_phi)] * exp(phi^)
		//注意，李代数上加法可交换，没有左右之分。扰动发生在李群上
		Vector_Add(phi,  delta_phi, 3,phi_1);
		//Exp_V3(phi_1,R);
		Rotation_Vector_3_2_Matrix(phi_1, R);
		Disp(R, 3, 3, "exp(phi + delta_phi)");		//此时，是最 exp(phi + delta_phi)

		//算exp(Jl(phi) * delta_phi)
		Get_Jl_3(phi, Jl);
		//Disp(Jl, 3, 3, "Jl");
		Matrix_Multiply_3x1(Jl, delta_phi, Temp);
		//Disp(Temp, 3, 1, "Jl*phi");
		//Exp_V3(Temp, Rl);
		Rotation_Vector_3_2_Matrix(Temp, Rl);
		
		//算exp(phi)
		//Exp_V3(phi, R);
		Rotation_Vector_3_2_Matrix(phi, R);
		//算exp(Jl(-1)*delta_phi) * exp(phi)
		Matrix_Multiply_3x3(Rl, R, R);
		Disp(R, 3, 3, "exp(Jl(phi)*delta_phi) * exp(phi)");
	}


	{//右扰动验证
		_T phi[3] = { 1,2,3 },
			delta_phi[4] = { 4, 3, 1, 0.05 },
			phi_1[3], R[3 * 3], Rr[3 * 3], Temp[3 * 3];
		_T Jr[3 * 3];

		//首先讲delta 化为3维向量
		Normalize(delta_phi, 3, delta_phi);
		Rotation_Vector_4_2_3<_T>(delta_phi, delta_phi);

		//蜀商的公式，exp(delta_phi + phi) = exp[(J(phi)*delta_phi)] * exp(phi^)
		//注意，李代数上加法可交换，没有左右之分。扰动发生在李群上
		Vector_Add(phi, delta_phi, 3, phi_1);
		//Exp_V3(phi_1,R);
		Rotation_Vector_3_2_Matrix(phi_1, R);
		Disp(R, 3, 3, "exp(phi + delta_phi)");		//此时，是最 exp(phi + delta_phi)

		//算exp(Jr(phi) * delta_phi)
		Get_Jr_3(phi, Jr);
		//Disp(Jr, 3, 3, "Jr");
		Matrix_Multiply_3x1(Jr, delta_phi, Temp);
		Rotation_Vector_3_2_Matrix(Temp, Rr);
		//求exp(phi)
		Rotation_Vector_3_2_Matrix(phi, R);
		
		//exp(phi) * exp(Jr(phi)*deltaphi)
		Matrix_Multiply_3x3(R, Rr, R);
		//Disp(R, 3, 3, "exp(phi) * exp(Jr(phi)*deltaphi)");
	}

	{//验证exp(delta_phi)exp(phi) = exp(phi + Jl(-1)(phi) * delta_phi)
		_T phi[3] = { 1,2,3 },
			delta_phi[4] = { 4, 3, 1, 0.05 },
			phi_1[3], R[3 * 3], Rl[3 * 3], Temp[3 * 3];
		_T Jl_Inv[3 * 3];

		//首先讲delta 化为3维向量
		Normalize(delta_phi, 3, delta_phi);
		Rotation_Vector_4_2_3<_T>(delta_phi, delta_phi);

		//蜀商的公式，exp(delta_phi + phi) = exp[(J(phi)*delta_phi)] * exp(phi^)
		//注意，李代数上加法可交换，没有左右之分。扰动发生在李群上
		Rotation_Vector_3_2_Matrix(delta_phi, Rl);
		Rotation_Vector_3_2_Matrix(phi, R);
		Matrix_Multiply_3x3(Rl, R, R);
		Disp(R, 3, 3, "R");

		//求J(-1)(phi)
		Get_Jl_Inv_3(phi, Jl_Inv);
		//Disp(Jl_Inv, 3, 3, "Jl_inv");
		//Jl(-1)(phi)*delta_phi
		Matrix_Multiply_3x1(Jl_Inv, delta_phi, phi_1);
		//Disp(phi_1, 3, 1, "Jl(-1)(phi)*delta_phi");
		Vector_Add(phi, phi_1,3, Temp);
		//求 exp(phi + J(-1)(phi)*delta_phi)
		Exp_V3(Temp, R);
		Disp(R, 3, 3, "R");
	}
	return;
}

void J_Test_1()
{//这个试验要揭示为什么李代数（旋转向量）不能直接加
	typedef double _T;
	_T P[3] = { 0,1,0 },							//车鼻子上一点
		V1[4] = { 0,0,1,fAngle_2_Radian<_T>(70) },		//先左转70度
		V2[4] = { 1, 0, 0, fAngle_2_Radian<_T>(10) };		//抬头10度

	//真实的物理世界
	_T Rl[3 * 3], R[3 * 3],P1[3];
	Rotation_Vector_4_2_Matrix(V1, R);
	Rotation_Vector_4_2_Matrix(V2, Rl);
	Matrix_Multiply_3x3(Rl, R, R);
	Matrix_Multiply_3x1(R, P, P1);
	Disp(P1, 3, 1, "Actural");

	//错误的做法，直接在李代数上加
	_T Temp[4], V3[4];
	Rotation_Vector_4_2_3(V1, V3);
	Rotation_Vector_4_2_3(V2, Temp);
	Vector_Add(Temp, V3, 3, V3);
	Exp_V3(V3, R);
	Matrix_Multiply_3x1(R, P, P1);
	Disp(P1, 3, 1, "Wrong");

	//近似
	_T Jl_Inv[3 * 3];
	Get_Jl_Inv_4(V1, Jl_Inv);
	Rotation_Vector_4_2_3(V2, Temp);

	Matrix_Multiply_3x1(Jl_Inv, Temp, Temp);
	Rotation_Vector_4_2_3(V1, V3);
	Vector_Add(V3, Temp,3, V3);
	Exp_V3(V3, R);
	Matrix_Multiply_3x1(R, P, P1);
	Disp(P1, 3, 1, "Correct");
	return;
}
void Deriv_Test_1()
{//求导验算
	typedef double _T;
	{//dRP/dphi
		_T V[4] = { 1,2,3,PI / 5 }, R[3 * 3];
		Normalize(V, 3, V);
		Rotation_Vector_4_2_Matrix(V, R);
		
		_T P[3] = { 2,3,4 }, Rp[3], dRP_dphi[3 * 3];
		Matrix_Multiply_3x1(R, P, Rp);

		Get_dRP_dphi(Rp, dRP_dphi);
		Disp(dRP_dphi, 3, 3, "dRP/dphi");		
	}
}


void QR_Real_Test_1()
{//
	typedef double _T;
	_T A[9] = { 1, 2, 3, 4, 5, 6, 7, 8, 9 }, R[9];

	int iResult = QR_Decompose_Real<_T>(A, 3, R);
	for (int i = 0; i < 3; i++)
	{
		_T fEigen = R[i * 3 + i];
		_T Eigen_Vector[3 * 3];
		int iCount = 0;
		Get_Eigen_Vector_by_Value_Real<_T>(A, 3, fEigen, Eigen_Vector,&iCount);
		Disp(Eigen_Vector, iCount, 3, "Eigen_Vector");
	}
	return;
}
void QR_Real_Test_2()
{
	typedef double _T;
	const int n = 10;
	_T  A[] = {-1293.8,      187.925,    347.56,     248.025,  -1508.625,   -443.58,   -1638.4965,  1485.87314,   -0.685,    1894.9116,
	187.925,     -7.9565,   -38.17445,  -31.565,    133.95,      70.173,     171.7936,   -154.88815,   26.157,    -190.7216,
	347.56,     -38.17445,   30.943,    -35.3135,   -210.91,     141.199,     -43.5535,     54.6465,   171.107,     40.265,
	248.025,    -31.565,    -35.3135,   -32.415,    124.385,     91.2525,    183.3117,   -159.251,     41.8495,   -216.5475,
  -1508.625,    133.95,     -210.91,    124.385,   1455.9,       -127.71,     654.5275,   -673.505,   -912.725,   -620.835,
   -443.58,      70.173,     141.199,    91.2525,   -127.71,     -148.23,    -649.8716,    593.9195,    34.375,     752.018,
  -1638.4965,   171.7936,   -43.5535,   183.3117,   654.5275,   -649.8716,   -24.457,       -41.7235, -730.12775,  136.6325,
   1485.87314, -154.88815,   54.6465,  -159.251,    -673.505,    593.9195,    -41.7235,     92.188,    683.054,   -253.53,
	 -0.685,     26.157,    171.107,     41.8495,   -912.725,     34.375,   -730.12775,    683.054,    250.92,     805.7865,
   1894.9116,   -190.7216,   40.265,   -216.5475,   -620.835,    752.018,     136.6325,   -253.53,     805.7865,   -282.13
	};

	_T R[n];
	//从这个例子看出，传统的QR分解几十实数有解，也会迭代太慢
	int iResult = QR_Decompose_Real(A, n, R);
	for (int i = 0; i < n; i++)
	{
		_T Eigen_Vector[n];
		_T fEigen = R[i];
		int iCount;
		Get_Eigen_Vector_by_Value_Real(A, n, fEigen, Eigen_Vector, &iCount);
		Disp(Eigen_Vector, iCount, n, "Eigen_Vector");
	}
	return;
}

void QR_Complex_Test_1()
{
	typedef double _T;
	_T A[] = { 0.0000, 1.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
		0.0000, 0.0000, 1.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
		0.0000, 0.0000, 0.0000, 1.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
		0.0000, 0.0000, 0.0000, 0.0000, 1.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
		0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 1.0000, 0.0000, 0.0000, 0.0000, 0.0000,
		0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 1.0000, 0.0000, 0.0000, 0.0000,
		0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 1.0000, 0.0000, 0.0000,
		0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 1.0000, 0.0000,
		0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 1.0000,
		0.1473, -0.0875, -0.1759, 0.1704, -0.5253, 1.0418, -1.6833, 2.7600, -2.1000, 1.7000  };
	int n = 10;
	Complex<_T> R[10], Q[10 * 10];
	int iReal_Count;
	int iResult = QR_Decompose_Complex(A, 10, R, &iReal_Count, Q);

	for (int i = 0; i < n; i++)
		Disp_Complex<_T>(R[i]);
	
	return;
}
int main_Slam_Test_2()
{
	bInit_Env_CPU(15000000);
	//Rotation_Test_1();
	//T_Inv_Test();
	//Draw_Pose_Test();
	//Exp_Test_1();
	//J_Test_1();
	//Chess_Board_Detect_Main();
	//Deriv_Test_1();
	//QR_Real_Test_2();
	QR_Complex_Test_1();
	Free_Env_CPU();
#ifdef WIN32
	_CrtDumpMemoryLeaks();
#endif
	return 0;
}