#include "slam.h"

static void E_Test_1()
{//35点对估计E
	typedef double _T;
	_T A[][2] = {
		-0.09182510, -0.02740728,
-0.08761927, -0.02519430,
0.09213333, -0.02308262,
0.09749971, -0.02258261,
-0.08129062, -0.00879808,
0.13106778, -0.00098159,
0.16098021, 0.00250948,
-0.09757352, 0.00503549,
0.14530507, 0.00761494,
0.09857728, 0.00814312,
0.06985272, 0.01174567,
-0.10424798, 0.01265123,
-0.02462738, 0.01558213,
-0.00950911, 0.01573378,
0.16002843, 0.01700515,
0.07011022, 0.01858117,
0.10259226, 0.04260552,
0.04726539, 0.05604459,
0.01003314, 0.06946079,
-0.05945094, 0.07189213,
0.17886913, -0.05598897,
0.14971074, -0.01989804,
-0.08225460, -0.01770928,
0.09877884, -0.00846430,
0.11130768, -0.00840081,
0.11515167, -0.00050014,
0.01795247, 0.01843837,
0.12437593, 0.02074700,
0.07700882, 0.02791510,
0.08111150, 0.02967108,
0.09300036, 0.03527236,
0.07288723, 0.06279132,
0.05312944, 0.08085094,
-0.08553912, -0.04613675,
-0.08703543, -0.01454860,
0.16512555, -0.00142763,
0.13580646, 0.00172531,
0.16006631, 0.00982405,
0.09722903, 0.02985964,
0.11713843, 0.04959117,
-0.09901395, -0.05313541,
-0.07525602, -0.02743964,
0.16047275, -0.00733985,
-0.02037770, 0.02243260,
0.06981742, 0.08364803,
0.07886662, 0.07345331,
0.16760470, -0.07789481,
0.11425643, 0.04072527};
	
	_T B[][2] = {
			-0.08938950, -0.04679362,
-0.08527513, -0.04481430,
0.09237506, -0.01970226,
0.09748751, -0.01993329,
-0.08064144, -0.02949643,
0.13081363, 0.00060231,
0.16110212, 0.00719443,
-0.09828678, -0.00742442,
0.14524268, 0.00964363,
0.09814857, 0.00919917,
0.06915431, 0.01233628,
-0.10591899, -0.00062450,
-0.02647541, 0.01834036,
-0.01106944, 0.01806002,
0.15983000, 0.02210762,
0.06949040, 0.01889808,
0.10110219, 0.04536230,
0.04468023, 0.06197237,
0.00526494, 0.07164607,
-0.06619190, 0.06925484,
0.17859043, -0.05307224,
0.14951219, -0.01556544,
-0.08111998, -0.03820267,
0.10217904, -0.00602620,
0.11150852, -0.00625814,
0.11512687, 0.00106496,
0.01654689, 0.02143814,
0.12407864, 0.02115393,
0.07568007, 0.02908716,
0.08027240, 0.03043711,
0.09194091, 0.03601122,
0.07069383, 0.06416433,
0.04858549, 0.08723086,
-0.08142594, -0.06309575,
-0.08570673, -0.03549355,
0.16504072, 0.00406775,
0.13559248, 0.00258636,
0.16027127, 0.01384947,
0.09649882, 0.03057039,
0.11634257, 0.05944812,
-0.09399816, -0.06886531,
-0.07246798, -0.04550490,
0.16047652, -0.00210445,
-0.02236267, 0.02892585,
0.06564952, 0.09160234,
0.07544430, 0.08139120,
0.16700771, -0.07748599,
0.11503830, 0.04825788};
	
	int iInlier,iCount = sizeof(A) / (sizeof(_T)*2);
	_T T[3 * 4], E[3 * 3], eps = (4.f / 768) * (4.f / 768);

	int iResult;
	if (!(iResult = Estimate_E_T_nPoint<_T>(A, B, iCount, &iInlier, T, E,eps,1)))
		return;

	Compute_Squared_Sampson_Error<_T>(E, A, B, iCount, &iInlier, eps);
	iInlier = Get_Inlier_Count(A, B, iInlier, T);
	//Test_E(E1, A, B, iInlier);
	//Disp(T, 3, 4, "T");
	//Disp(E, 3, 3, "E");
	printf("%d inlier found\n", iInlier);
	return;
}

void Sift_Test_1()
{
	Point_2D<float>* pKey_Point;
	unsigned char (*pFeature)[128];
	int iCount;
	Get_Sift_Feature("C:\\tmp\\dev_env\\sample\\01.bmp", &pKey_Point, &pFeature, &iCount, 0);
	Free(pKey_Point);
	Free(pFeature);
}

void Sift_Test_2()
{
	typedef double _T;
	Pose_Base<_T> oPose_Base;
	if (!Pose_Base_Init(&oPose_Base))
		return;
	float (*pKey_Point)[2];
	unsigned char (*pFeature)[128];
	int iResult, iCount, iID;

	unsigned long long tStart = iGet_Tick_Count();

	{
		//获得Keypoint
		iResult = Get_Sift_Feature_1("c:\\tmp\\dev_env\\Sample\\01.bmp", &pKey_Point, &pFeature, &iCount, 0);
		//加入Base
		Pose_Base_Add_Feature<_T>(&oPose_Base, pKey_Point, pFeature, iCount, &iID);
		Free(pKey_Point);
		Free(pFeature);
	}

	{
		iResult = Get_Sift_Feature_1("c:\\tmp\\dev_env\\Sample\\02.bmp", &pKey_Point, &pFeature, &iCount, 0);
		Pose_Base_Add_Feature<_T>(&oPose_Base, pKey_Point, pFeature, iCount, &iID);
		Free(pKey_Point);
		Free(pFeature);
	}

	{
		unsigned char (*pFeature_0)[128], (*pFeature_1)[128];
		float (*pKey_Point_0)[2], (*pKey_Point_1)[2];
		int iCount_0, iCount_1;
		pFeature_0 = Pose_Base_Get_Feature(&oPose_Base, 0, &iCount_0);
		pFeature_1 = Pose_Base_Get_Feature(&oPose_Base, 1, &iCount_1);
		pKey_Point_0 = Pose_Base_Get_Key_Point(&oPose_Base, 0, NULL);
		pKey_Point_1 = Pose_Base_Get_Key_Point(&oPose_Base, 1, NULL);

		unsigned short (*pMatch_Pair)[2];
		int iMatch_Count;
		Sift_Match_2_Image(pKey_Point_0, pFeature_0, iCount_0,
			pKey_Point_1, pFeature_1, iCount_1, &pMatch_Pair, &iMatch_Count);

		//Sift_Match_2_Image(pFeature_0, iCount_0, pFeature_1, iCount_1,&pMatch_Pair,&iMatch_Count);
		Pose_Base_Add_Match(&oPose_Base, 0, 1, pMatch_Pair, iMatch_Count);

		_T(*pPoint_A)[2], (*pPoint_B)[2];
		Pose_Base_Get_Match_Key_Point(&oPose_Base, 0, 1, &pPoint_A, &pPoint_B, &iCount);

		_T K[] = { 1520.4,1525.9,302.32,246.87 },
			//K[]={ 1523.15 ,1523.15 ,9,302.32,246.87 },
			D[5] = { 0 }, T[3 * 4], E[3 * 3];

		uv_2_Norm(pPoint_A, iCount, K, 4, pPoint_A);
		uv_2_Norm(pPoint_B, iCount, K, 4, pPoint_B);
		//Disp((_T*)pPoint_A, iMatch_Count, 2, "A");
		//Disp((_T*)pPoint_A, iMatch_Count, 2, "A");
		int iInlier;
		_T fError;
		Ransac_E<_T>(pPoint_A, pPoint_B, iCount, T, E, &iInlier,&fError, 1, 0);
		//int iInlier = Get_Inlier_Count(pPoint_A, pPoint_B, iCount, T, 1, 1, &fError);
		printf("Inlier:%d\n", iInlier);
		//Disp((_T*)pPoint_A, iCount,2, "A");
		//Disp((_T*)pPoint_B, iCount, 2, "B");
	}
	{
		float(*pPoint_A)[2], (*pPoint_B)[2];
		int iCount;
		pPoint_A = Pose_Base_Get_Key_Point(&oPose_Base, 0, &iCount);
		pPoint_B = Pose_Base_Get_Key_Point(&oPose_Base, 1, &iCount);
		//Disp((float*)pPoint_B, iCount, 2,"A");
		//Ransac_E<_T>(pPoint_A, pPoint_B, NULL, NULL, 1, 0);
	}
	/*{
		const int iBlock_Width = 32;
		unsigned char Source[iBlock_Width * iBlock_Width], Dest[iBlock_Width * iBlock_Width];
		for (int y = 0; y < iBlock_Width; y++)
			for (int x = 0; x < iBlock_Width; x++)
				Source[y * iBlock_Width + x] = x;
		Disp(Source, iBlock_Width, iBlock_Width, "Source");
		Transpose_32x32(Source, iBlock_Width, Dest, iBlock_Width);
		Disp(Dest, iBlock_Width, iBlock_Width, "Dest");
	}*/

	Pose_Bsase_Free(&oPose_Base);
	return;
}

void Two_Image_Pose_Test()
{
	typedef float _T;
	_T T[3*4],K[] = { 1520.4f,1525.9f,302.32f,246.87f }, D[5] = { 0 };
	_T (*pA)[2], (*pB)[2],fError;
	int iResult, iCount;
	
	{
		iResult = Get_2_Image_Pose_by_File<_T>(
			"c:\\tmp\\dev_env\\sample\\01.bmp",
			"c:\\tmp\\dev_env\\sample\\02.bmp",
			K, T, &pA, &pB, &iCount, &fError);
		if (!iResult)
			return;
		Free(pA);
	}

	{
		iResult = Get_2_Image_Pose_by_File<_T>("c:\\tmp\\dev_env\\sample\\02.bmp",
			"c:\\tmp\\dev_env\\sample\\03.bmp", K, T, &pA, &pB, &iCount, &fError,0);
		if (!iResult)
			return;

		//验算
		iResult = Get_Inlier_Count(pA, pB, iCount, T);
		Free(pA);
	}
	return;
}

void Pose_Base_Test()
{

}

int Sift_Main()
{
	bInit_Env_CPU(40000000);
	//Pose_Base_Test();
	Two_Image_Pose_Test();
	//E_Test_1();
	//Sift_Test_2();
	//Gauss_Filter_Test_1();
	Free_Env_CPU();
#ifdef WIN32
	_CrtDumpMemoryLeaks();
#endif
	return 0;
}