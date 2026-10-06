//五点法估计E矩阵
#include "Slam.h"
using namespace std;

void Poly_Test_1()
{//多项式的表示
    typedef double _T;
    {
        Polynormial<_T> oA, oB, oA_B;
        Init_Poly(&oA, 2, 100);
        Init_Poly(&oB, 2, 100);
        Insert_Term(&oA, oGet_Term_7<_T>(1.2, 3, 3));    //A = 1.2* x^2
        Insert_Term(&oA, oGet_Term_7<_T>(3.4, 3, 2));    //+ 3.4 * y

        Insert_Term(&oB, oGet_Term_7<_T>(5.6, 2, 3));    //A = 1.2* x^2
        Insert_Term(&oB, oGet_Term_7<_T>(7.8, 0, 0));    //A = 1.2* x^2
        //Disp_Poly(oA, "A");
        //Disp_Poly(oB, "B");
        oA_B = oA + oB;
        //Disp_Poly(oA_B, "C");

        Free_Poly(&oA);
        Free_Poly(&oB);
        Free_Poly(&oA_B);
    }

    {//多项累加
        Polynormial<_T> oA, oB;
        Init_Poly(&oA, 2, 100);
        Init_Poly(&oB, 2, 100);
        Insert_Term(&oA, oGet_Term_7<_T>(1.2, 3, 1, 2, 1, 1, 3));   //a = x^3 + y^2 + z^3
        Insert_Term(&oA, oGet_Term_7<_T>(1.2, 3, 3, 2, 2, 1, 1));   
        //Disp_Poly(oA, "A");
        Free_Poly(&oA);
        Free_Poly(&oB);
    }

    {
        Polynormial<_T> oA, oB, oC;
        Init_Poly(&oA, 2, 100);
        Init_Poly(&oB, 2, 100);
        Insert_Term(&oA, oGet_Term_7<_T>(1));
        Insert_Term(&oA, oGet_Term_7<_T>(1, 1,1,3,2));

        Insert_Term(&oB, oGet_Term_7<_T>(1 ));
        Insert_Term(&oB, oGet_Term_7<_T>(1, 2,1,3,1));
        Disp_Poly(oA, "A");
        Disp_Poly(oB, "B");
        oC = oA * oB;
        Disp_Poly(oC, "C");
        Free_Poly(&oA);
        Free_Poly(&oB);
        Free_Poly(&oC);
    }
    return;
}

void Poly_Test_2()
{
    typedef double _T;
    _T E1[3 * 3], E2[3 * 3], E3[3 * 3], E4[3 * 3];
    for (int y = 0; y < 3; y++)
        for (int x = 0; x < 3; x++)
            E1[y * 3 + x] = (_T)1 * 10.f + y + (_T)x / 10,
            E2[y * 3 + x] = (_T)2 * 10.f + y + (_T)x / 10,
            E3[y * 3 + x] = (_T)3 * 10.f + y + (_T)x / 10,
            E4[y * 3 + x] = (_T)4 * 10.f + y + (_T)x / 10;

    Polynormial<_T> E[3 * 3];

    for(int y=0;y<3;y++)
    {
        for (int x = 0; x < 3; x++)
        {
            int iPos = y * 3 + x;
            Init_Poly(&E[iPos], 4, 4);
            Poly_Term_7<_T> oTerm = oGet_Term_7(E1[iPos], 1, 1);
            Insert_Term(&E[iPos], oTerm);

            oTerm = oGet_Term_7(E2[iPos], 2, 1);
            Insert_Term(&E[iPos], oTerm);
            
            oTerm = oGet_Term_7(E3[iPos], 3, 1);
            Insert_Term(&E[iPos], oTerm);

            oTerm = oGet_Term_7(E4[iPos]);
            Insert_Term(&E[iPos], oTerm);
            Disp_Poly(E[iPos], "Poly");
        }
    }
    
    //第一堆： E00*(E11*E22 - E12*E21)
    Polynormial<_T> E11_E22 = E[1 * 3 + 1] * E[2 * 3 + 2];
    Disp_Poly(E11_E22, "E11*E22");

    Polynormial<_T> E12_E21 = E[1 * 3 + 2] * E[2 * 3 + 1];
    Disp_Poly(E12_E21, "E12*E21");

    Polynormial<_T>E11_E22_E12_E21 = E11_E22 - E12_E21;
    Disp_Poly(E11_E22_E12_E21, "E11*E22 - E12*E21");

    Polynormial<_T> E00_E11_E22_E12_E21 = E[0 * 3 + 0] * E11_E22_E12_E21;
    Disp_Poly(E00_E11_E22_E12_E21, "E00*(E11*E22 - E12*E21)");

    //第二坨
    Polynormial<_T> E10_E22 = E[1 * 3 + 0] * E[2 * 3 + 2];
    Disp_Poly(E10_E22, "E10*E22");

    Polynormial<_T>E12_E20 = E[1 * 3 + 2] * E[2 * 3 + 0];
    Disp_Poly(E12_E20, "E12*E20");

    Polynormial<_T> E10_E22_E12_E20 = E10_E22 - E12_E20;
    Disp_Poly(E10_E22_E12_E20, "E10*E22-E12*E20");

    Polynormial<_T> E01_E10_E22_E12_E20 = E[0 * 3 + 1] * E10_E22_E12_E20;
    Polynormial<_T> oNeg = E01_E10_E22_E12_E20 * ((_T)-1);
    Disp_Poly(E01_E10_E22_E12_E20, "E01*(E10*E22-E12*E20)");
    Disp_Poly(oNeg, "-E01*(E10*E22-E12*E20)");

    //第三坨
    Polynormial<_T>E10_E21 = E[1 * 3 + 0] * E[2 * 3 + 1];
    Disp_Poly(E10_E21, "E10*E21");

    Polynormial<_T>  E11_E20 = E[1 * 3 + 1] * E[2 * 3 + 0];
    Disp_Poly(E11_E20, "E11*E20");

    Polynormial<_T> E10_E21_E11_E20 = E10_E21 - E11_E20;
    Disp_Poly(E10_E21_E11_E20, "E10*E21-E11*E20");

    Disp_Poly(E[0 * 3 + 2], "E02");

    Polynormial<_T>E02_E10_E21_E11_E20 = E[0 * 3 + 2] * E10_E21_E11_E20;
    Disp_Poly(E02_E10_E21_E11_E20, "E02*(E10*E21-E11*E20)");

    //三项相加
    Polynormial<_T> Sum_1 = E00_E11_E22_E12_E21 + E02_E10_E21_E11_E20;
    Polynormial<_T> Sum_2 = Sum_1 - E01_E10_E22_E12_E20;
    //Disp_Poly(Sum_1, "Det");
    //Disp_Poly(Sum_2, "E01");

    //第二堆，算 E01 * 
    for (int i = 0; i < 9; i++)
        Free_Poly(&E[i]);
    return;
}

void Inv_Test_1()
{//改造求逆过程，分离出基本行变换算法
    typedef double _T;
    _T A[3 * 3] = { 1,2,3,3,2,1,9,7,8 }, A_Inv[3 * 3];
    int iResult = Get_Inv_Matrix_Row_Op(A, A_Inv, 3);
    Disp(A_Inv, 3, 3, "Inv");
    _T fError;
    Test_Inv_Matrix(A, A_Inv, 3,&fError);

    return;
}

void E_Test_5()
{
    typedef double _T;
    //_T Match_Point[2][5][2] = {
    //    {
    //        { 0.244484,  0.146945},  // 点 1
    //        {-0.157251,  0.235690},  // 点 2
    //        { 0.108349, -0.271299},  // 点 3
    //        {-0.099707, -0.066426},  // 点 4
    //        { 0.177726,  0.118632}   // 点 5
    //    },
    //    {
    //        { 0.496925, -0.043578},  // 点 1
    //        { 0.103552,  0.002157},  // 点 2
    //        { 0.465784, -0.404004},  // 点 3
    //        { 0.162066, -0.223655},  // 点 4
    //        { 0.426421, -0.083895}   // 点 5
    //    }
    //};

    //// 满血版 5 参数畸变系数 (直接拷贝进代码)
    ////_T Distort[5] = { -0.250000,0.120000,-0.050000,0.001500,-0.002500 };

    //_T Match_Point[2][5][2] = {
    //{
    //    {  0.099623,  0.049846},  // 点 1
    //    { -0.055388,  0.121782},  // 点 2
    //    {  0.054335, -0.090589},  // 点 3
    //    { -0.031244, -0.024990},  // 点 4
    //    {  0.073544,  0.010517},  // 点 5
    //},
    //{
    //    {  0.417659, -0.074158},  // 点 1
    //    {  0.278874, -0.005762},  // 点 2
    //    {  0.361586, -0.212198},  // 点 3
    //    {  0.318330, -0.153758},  // 点 4
    //    {  0.398155, -0.114513},  // 点 5
    //}
    //};
    //_T Distort[5] = { -0.250000, 0.120000, -0.050000, 0.001500, -0.002500 };

    _T Match_Point[2][5][2] = {
        {
            { 0.049331,  0.097400},  // 点 1
            { 0.134389, -0.170512},  // 点 2
            { 0.088206,  0.084746},  // 点 3
            {-0.078250,  0.059111},  // 点 4
            {-0.048556, -0.005626},  // 点 5
        },
        {
            { 0.359063,  0.020132},  // 点 1
            { 0.449098, -0.196788},  // 点 2
            { 0.439679,  0.020969},  // 点 3
            { 0.180286,  0.002138},  // 点 4
            { 0.176803, -0.054028},  // 点 5
        }
    };
    _T Distort[5] = { -0.250000, 0.120000, -0.050000, 0.001500, -0.002500 };

    //先对点集去畸变
    for (int i = 0; i < 5; i++)
    {
        if (!Get_Pn_by_Pd(Match_Point[0][i], Distort, Match_Point[0][i]))
            printf("error");
        if (!Get_Pn_by_Pd(Match_Point[1][i], Distort, Match_Point[1][i]))
            printf("Error");

        //Disp(Match_Point[0][i], 1, 2,"A");
        //Disp(Match_Point[1][i], 1, 2, "B");
    }

    _T R[3 * 3], t[3],E[3*3];
    int iResult = Estimate_E_Rt_5Point<_T>(Match_Point[0], Match_Point[1],R,t,E);
    Disp(E, 3, 3, "E");
    Disp(R, 3, 3, "R");
    Disp(t, 1, 3, "t");

    _T T[4 * 4];
    Gen_Pose_By_R_t(R, t, T);
    //Test_E(E, Match_Point[0], Match_Point[1], 5);

    unsigned long long tStart = iGet_Tick_Count();
    _T fError;
    //for(int i=0;i<1000000;i++)
        fError = Test_Triangulate(T, Match_Point[0], Match_Point[1], 5);

    printf("%e %lld\n", fError, iGet_Tick_Count() - tStart);;
    return;
}

void E_Test_6()
{
    typedef double _T;;
    _T R[9] = { 1.000000, 0.000000, 0.000000, 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 1.000000 };
    _T t[3] = { 0.000000, 0.000000, 0.120000 };

    _T Match_Point[][1][2] = { -0.018390, 0.034617,-0.029096, 0.019554 };
        
    _T T[4*4],fError;
    Gen_Pose_By_R_t(R, t, T);
    unsigned long long tStart = iGet_Tick_Count();
    for(int i=0;i<1000000;i++)
    fError = Test_Triangulate(T, Match_Point[0], Match_Point[1], 1);
    printf("Error:%f %lld\n", fError, iGet_Tick_Count() - tStart);;
    return;
}

int main()
{
	bInit_Env_CPU(15000000,128);
    E_Test_6();
	Free_Env_CPU();
#ifdef WIN32
	_CrtDumpMemoryLeaks();
#endif
	return 0;
}