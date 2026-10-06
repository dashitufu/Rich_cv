#include "Slam.h"



void E_Test_1()
{
    typedef double _T;
    _T Match_Point[2][8][2] = {
        {
            { -0.10145773,  0.03172393},  // 视1 点1
            {  0.01233337, -0.12288783},  // 视1 点2
            {  0.13829576,  0.14016054},  // 视1 点3
            {  0.08234948, -0.06507584},  // 视1 点4
            { -0.14983571, -0.07444965},  // 视1 点5
            { -0.27236866,  0.08109024},  // 视1 点6
            {  0.01061068,  0.12990613},  // 视1 点7
            { -0.14262515, -0.11262459},  // 视1 点8
        },
        {
            {  0.14176422,  0.03773067},  // 视2 点1
            {  0.16822967, -0.12437105},  // 视2 点2
            {  0.49335097,  0.15835162},  // 视2 点3
            {  0.51871158, -0.04721235},  // 视2 点4
            {  0.17490612, -0.06799633},  // 视2 点5
            {  0.12715178,  0.10159161},  // 视2 点6
            {  0.20945628,  0.13520051},  // 视2 点7
            {  0.15399902, -0.10898468},  // 视2 点8
        }
    };
    
    _T E[3 * 3], T[4 * 4];
    int iResult;
    unsigned long long tStart = iGet_Tick_Count();
    for(int i=0;i<100000;i++)
    iResult = Estimate_E_8Point<_T>(Match_Point[0], Match_Point[1], T, E);
    printf("%lld\n", iGet_Tick_Count() - tStart);
    //Disp(T, 3, 4, "T");
    //Disp(E, 3, 3, "E");
    printf("Error: %e\n", Test_Triangulate(T, Match_Point[0], Match_Point[1], 8));

    return;
}
int main()
{
    bInit_Env_CPU(15000000, 128);
    E_Test_1();
    Free_Env_CPU();
#ifdef WIN32
    _CrtDumpMemoryLeaks();
#endif
    return 0;
}