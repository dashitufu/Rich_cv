#include "Slam.h"
#include "Chess_Board_Detect.h"

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

void F_Test_1()
{//试一下七点发
    typedef double _T;
    _T Match_Point[2][7][2] = {
        {
            {    744.2859,   173.3696},  // 视1 点1
            {    343.8562,   532.5118},  // 视1 点2
            {    914.9771,   511.5083},  // 视1 点3
            {    439.8113,   243.5585},  // 视1 点4
            {    658.3596,   627.7341},  // 视1 点5
            {    863.9284,   246.4212},  // 视1 点6
            {    521.7593,   354.4078},  // 视1 点7
        },
        {
            {   1038.6567,   158.2253},  // 视2 点1
            {    593.1858,   504.5274},  // 视2 点2
            {   1158.6475,   506.5366},  // 视2 点3
            {    655.8702,   230.2871},  // 视2 点4
            {    938.0605,   611.7322},  // 视2 点5
            {   1092.8117,   227.5149},  // 视2 点6
            {    751.5577,   338.7041},  // 视2 点7
        }
    };
    _T F[3][3 * 3];
    int iCount;
    int iResult = Estimate_F_7Point<_T>(Match_Point[0], Match_Point[1], F,&iCount);
    return;
}

void F_Test_2()
{
    typedef double _T;
    _T Match_Point[2][8][2] = {
        {
            {   744.2859,   173.3696},  // 视1 点1
            {   343.8562,   532.5118},  // 视1 点2
            {   914.9771,   511.5083},  // 视1 点3
            {   439.8113,   243.5585},  // 视1 点4
            {   658.3596,   627.7341},  // 视1 点5
            {   863.9284,   246.4212},  // 视1 点6
            {   521.7593,   354.4078},  // 视1 点7
            {   201.6848,   114.5674},  // 视1 点8   <-- 新增
        },
        {
            {  1038.6567,   158.2253},  // 视2 点1
            {   593.1858,   504.5274},  // 视2 点2
            {  1158.6475,   506.5366},  // 视2 点3
            {   655.8702,   230.2871},  // 视2 点4
            {   938.0605,   611.7322},  // 视2 点5
            {  1092.8117,   227.5149},  // 视2 点6
            {   751.5577,   338.7041},  // 视2 点7
            {   509.7803,   116.3734},  // 视2 点8   <-- 新增
        }
    };
    _T F[3 * 3];
    unsigned long long tStart = iGet_Tick_Count();
    for(int i=0;i<100000;i++)
    Estimate_F_8Point(Match_Point[0], Match_Point[1], F);
    printf("%lld\n", iGet_Tick_Count() - tStart);
}


void H_Test_1()
{
    typedef double _T;
    _T Match_Point[2][4][2] = {
        {
            { 0.493026, 0.387309 },   // 点 1
            { 0.782311, 0.351190 },   // 点 2
            { 0.472040, 0.557427 },   // 点 3
            { 0.792034, 0.559527 },   // 点 4
        },
        {
            { 0.513000, 0.403000 },   // 点 1
            { 0.841775, 0.377884 },   // 点 2
            { 0.498118, 0.588222 },   // 点 3
            { 0.867903, 0.613124 },   // 点 4
        }
    };
    _T Distort[5] = { -0.100000, 0.020000, 0.000000, 0.000000, 0.000000 };
    for(int i=0;i<4;i++)
    {
        Get_Pn_by_Pd(Match_Point[0][i], Distort, Match_Point[0][i]);
        Get_Pn_by_Pd(Match_Point[0][i], Distort, Match_Point[1][i]);
    }

    //Disp((_T*)Match_Point[0], 4, 2, "A");
    //Disp((_T*)Match_Point[1], 4, 2, "B");

    _T H[3 * 3],T[3*4];
    int iResult;
    unsigned long long tStart = iGet_Tick_Count();
    for(int i=0;i<100000;i++)
    iResult = Estimate_H_4Point(Match_Point[0], Match_Point[1], T, H);
    printf("%lld\n", iGet_Tick_Count() - tStart);

    return;
}

int main()
{
    bInit_Env_CPU(15000000, 128);
    H_Test_1();
    //Chess_Board_Detect_Main();
    //Free_Env_CPU();
    Free_Env_CPU();
#ifdef WIN32
    _CrtDumpMemoryLeaks();
#endif
    return 0;
}