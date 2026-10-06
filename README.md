# Rich_cv
My Computer Vision Project

1，8点法估计E矩阵

    int iInlier, iResult = Ransac_E(Match_Point[0], Match_Point[1], 200, T, E, 0);
    iInlier = Get_Inlier_Count(Match_Point[0], Match_Point[1], 200, T, 1, 2, &fError);

2，5点法估计E矩阵

    int iInlier, iResult = Ransac_E(Match_Point[0], Match_Point[1], 200, T, E, 1);
    iInlier = Get_Inlier_Count(Match_Point[0], Match_Point[1], 200, T, 1, 2, &fError);
   
