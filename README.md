# Rich_cv
My Computer Vision Project

1，8点法估计E矩阵

    int iInlier, iResult = Ransac_E(Match_Point[0], Match_Point[1], 200, T, E, 0);
    iInlier = Get_Inlier_Count(Match_Point[0], Match_Point[1], 200, T, 1, 2, &fError);

2，5点法估计E矩阵

    int iInlier, iResult = Ransac_E(Match_Point[0], Match_Point[1], 200, T, E, 1);
    iInlier = Get_Inlier_Count(Match_Point[0], Match_Point[1], 200, T, 1, 2, &fError);
   
3，4点法估计H矩阵

    int iResult = Ransac_H<_T>(Match_Point[0], Match_Point[1], 200, T, H, eps);
    iInlier = Get_Inlier_Count(Match_Point[0], Match_Point[1], 200, T, 1, 2, &fError, eps);

4，两图估计位姿，内点，投影误差

    iResult = Get_2_Image_Pose_by_File<_T>(
	"c:\\tmp\\dev_env\\sample\\01.bmp",
	"c:\\tmp\\dev_env\\sample\\02.bmp",
	K, T, &pA, &pB, &iCount, &fError);
    
