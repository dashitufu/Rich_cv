#include "Image.h"
#include "Slam.h"
#include <algorithm>
using namespace std;

#define VL_SIFT_BILINEAR_ORIENTATIONS 1
typedef struct Sift_Keypoint {//关键点
	short ix, iy;		// x,y
	char is;			//第几层
	char o;				//Current octave，第几轮
	float x, y, s;		//经过泰勒修正后的小数值
	float sigma;
}Sift_Keypoint;

#define MAX_KEY_POINT_COUNT 1280
typedef struct Sift_Feature {
	short x_Keypoint, y_Keypoint;	//原来在该分辨率下Keypoint的位置
	float angle;		//主方向，这个源代码没有，但是觉得有必要保留下来，对进一步分析有用
	float x, y;			//这个位置的计算尚有意义，感觉对不上
	short o;			//在哪一轮，从这一版开始，octave指的是图像分辨率轮次
	short s;			//在某一轮的哪一层。理论上是指梯度层？
	unsigned char m_Desc_i[128];
	//float a11, a12, a21, a22;	//这几个有异议
}Sift_Feature_Byte;

typedef struct Sift_Feature_Level {
	Sift_Feature* m_pFeature;
	Sift_Keypoint* m_pKeypoint;
	int m_iFeature_Count, m_iKeypoint_Count;
}Sift_Feature_Level;

typedef struct Sift {
	int O;               //一共可以搞多少轮，好像没有必要留
	int S;               //每轮多少层，也是可以设置
	int o_min;           //可以-1，也可以更高一级。比如4K的情况下没有必要搞那么大
	int s_min;           /**< minimum level index. */
	int s_max;           /**< maximum level index. */
	int o_cur;           /**< current octave. */

	//分析结果，一共有6层多分辨核高斯滤波，5层DoG, 3层梯度
	//故此，按理说内存分配3层DoG,3层Octave即可
	//尝试以下内存安排
	//|  /   | G_0 | DoG_0 | G_1 | DoG_1 | G_2
	//|DoG_2 | G_3 | DoG_3 | G_4 | DoG_4 | G_5
	float* m_pDoG[3];		//current DoG data, DoG= difference of Gaussians 高斯差
	//注意，在m_pDoG与高斯图之间还要加上一行，用作最后一行的缓冲
	float* m_pGauss[3];		//装不同分辨率的高斯图。为避免二义性，将原来octave改名

	float* m_pGrad;			//重新安排内存试一下，m_pGrad与Octave DoG共用
	float* m_pTemp;			//对应 f->temp，重用当层的m_pDoG

	int grad_o;					//GSS gradient data octave.
	int m_iWidth, m_iHeight;	//源图的长宽

	//以下几个Sigma是关于高斯图的模糊程度，具体再代码中，后买你再总结
	float sigman;       //nominal image smoothing. 
	float sigma0;       //smoothing of pyramid base. 
	float sigmak;       //k-smoothing 
	float dsigma0;      //delta-smoothing. 

	//float* Gauss_Filter;  /**< current Gaussian filter */
	//float gaussFilterSigma;   /**< current Gaussian filter std */
	//int gaussFilterWidth;  /**< current Gaussian filter width */

	int octave_width;    // current octave width
	int octave_height;   // current octave height.

	int keys_res;        // size of the keys buffer. 

	float peak_thresh;  // peak threshold. 
	float edge_thresh;  // edge threshold. 
	float norm_thresh;  // norm threshold. 
	float magnif;       // magnification factor. 
	float windowSize;   // size of Gaussian window (in spatial bins) 

	int nkeys;           // Keypoint Count
	Sift_Keypoint* m_pKey_Point;	//Keypoint只是候选点，并不等于Feature

	int m_iMax_Feature_Count;
	int m_iMax_Keypoint_Count;		
	int m_iFeature_Count;			//最终的Feature Count
	//此处存放最终的整数形式Feature
	Sift_Feature* m_pFeature;		//最终的Feature所在，要的只是这一点点结果

	unsigned char* m_pBuffer_Start;	//每一轮内存开始之处
}Sift;

template void Pose_Bsase_Free(Pose_Base<double>* poBase);
template<typename _T>void Pose_Bsase_Free(Pose_Base<_T>* poBase)
{
	Free(poBase->m_pFeature_Buffer);
	Free(poBase->m_pMatch_Index);
	Free(poBase->m_pMatch_Pose);
	Free(poBase->m_pPose);
	*poBase = {};
}

template int Pose_Base_Init(Pose_Base<double>* poPose_Base,
	int iMax_Pose_Count, int iMax_Key_Point_Count,
	int iMax_Match_Pose_Count, int iMax_Match_Key_Point_Count);
template<typename _T>int Pose_Base_Init(Pose_Base<_T> *poPose_Base,
	int iMax_Pose_Count, int iMax_Key_Point_Count,
	int iMax_Match_Pose_Count, int iMax_Match_Key_Point_Count)
{//初始化一个大的位姿库
	Pose_Base<_T> oBase = {};	//一股脑儿初始化为0

	oBase.m_iMax_Key_Point_Count = iMax_Key_Point_Count;
	int iSize = oBase.m_iMax_Key_Point_Count * ((sizeof(float) * 2) + 128);
	oBase.m_pFeature_Buffer = (unsigned char(*)[128])pMalloc(iSize);
	oBase.m_pKey_Point_Buffer =(float(*)[2])(oBase.m_pFeature_Buffer + oBase.m_iMax_Key_Point_Count);
	
	oBase.m_iMax_Match_Index_Count = iMax_Match_Key_Point_Count;
	oBase.m_pMatch_Index = (unsigned short(*)[2])pMalloc(iMax_Match_Key_Point_Count * 2 * sizeof(unsigned short));
	
	oBase.m_iMax_Match_Pose_Count = iMax_Match_Pose_Count;;
	oBase.m_pMatch_Pose = (Match_Pose*)pMalloc(iMax_Match_Pose_Count * sizeof(Match_Pose));
	memset(oBase.m_pMatch_Pose, 0, iMax_Match_Pose_Count * sizeof(Match_Pose));

	oBase.m_iMax_Pose_Count = iMax_Pose_Count;
	oBase.m_pPose = (Pose_Bsse_Pose<_T> *)pMalloc(iMax_Pose_Count * sizeof(Pose_Bsse_Pose<_T>));
	memset(oBase.m_pPose, 0, iMax_Pose_Count * sizeof(Pose_Bsse_Pose<_T>));

	if (!oBase.m_pFeature_Buffer ||
		!oBase.m_pMatch_Index ||
		!oBase.m_pMatch_Pose ||
		!oBase.m_pPose)
	{
		printf("Fail to alloca memory in Init_Pose_Base\n");
		Free(oBase.m_pFeature_Buffer);
		Free(oBase.m_pMatch_Index);
		Free(oBase.m_pMatch_Pose);
		Free(oBase.m_pPose);
		return 0;
	}

	*poPose_Base = oBase;
	return 1;
}

static unsigned int iGet_Sift_Detect_Size(int iWidth, int iHeight, int o_min)
{//检测特征点阶段，预计算一共要用多少内存
//o_min: 必须要考虑进来，一旦不是-1开始，内存就不一样
	unsigned int iSize = 0, iSize_1 = 0;
	//Feature
	//iSize = (MAX_KEY_POINT_COUNT * 2 * (iWidth * iHeight / 307200))* sizeof(Sift_Feature);
	//iSize = MAX_KEY_POINT_COUNT * 2 * sizeof(Sift_Feature);
	//iSize = ALIGN_SIZE_128(iSize);
	////Keypoint
	//iSize += MAX_KEY_POINT_COUNT * sizeof(Sift_Keypoint);
	//iSize = ALIGN_SIZE_128(iSize);

	if (o_min < 0)
	{
		iWidth <<= -(o_min);
		iHeight <<= -(o_min);
	}
	else if (o_min > 0)
	{
		iWidth >>= o_min;
		iHeight >>= o_min;
	}
	iSize_1 = iWidth * iHeight * sizeof(float);
	iSize_1 = ALIGN_SIZE_128(iSize_1);

	//Octave
	iSize += iSize_1 * 3;

	//DoG，最小三层
	iSize += iSize_1 * 3 + iWidth * 3 * sizeof(float);

	//多加一行，用来算Grad最后一行
	//iSize += iWidth * 2 * sizeof(float);

	//grad 与Temp 共享
	//iSize += iSize_1 * 2;
	iSize = ALIGN_SIZE_128(iSize) + 128;	//留有余地128字节

	return iSize;
}

static void Init_Sift(Sift* poSift, int iOctave_Count, int iLayer_Count,
	int o_min, int iWidth, int iHeight, float fPeak_Thresh = 2.0 / 300.f)
{//此处应该一次性分配好所有的内存
	Sift oSift = {};
	//int iSize;

	oSift.m_iWidth = iWidth;
	oSift.m_iHeight = iHeight;

	//此处决定关键点选取阀值，值越大，掂越少
	//oSift.peak_thresh = 2.f / 300.f;
	oSift.peak_thresh = fPeak_Thresh;
	oSift.edge_thresh = 10.f;

	oSift.O = 4;
	oSift.S = 3;
	oSift.o_min = o_min;
	oSift.s_min = -1;
	oSift.s_max = iLayer_Count + 1;
	oSift.o_cur = oSift.o_min;

	oSift.sigman = 0.5;
	oSift.sigmak = (float)pow(2.0, 1.0 / iLayer_Count);
	oSift.sigma0 = 1.6f * oSift.sigmak;
	oSift.dsigma0 = (float)(oSift.sigma0 * sqrt(1.0f - 1.0f / (oSift.sigmak * oSift.sigmak)));

	oSift.octave_width = 0;
	oSift.octave_height = 0;

	oSift.nkeys = 0;

	oSift.norm_thresh = 0.0;
	oSift.magnif = 3.0;
	oSift.windowSize = 4 / 2;

	////安排内存
	//unsigned char* pCur = pBuffer;
	oSift.m_iMax_Feature_Count = oSift.m_iMax_Keypoint_Count= MAX_KEY_POINT_COUNT;	// *2 * (iWidth * iHeight / 307200);
	//iSize = oSift.m_iMax_Feature_Count * sizeof(Sift_Feature);
	//iSize = ALIGN_SIZE_128(iSize);

	//oSift.m_iFeature_Count = 0;
	//oSift.m_pFeature = (Sift_Feature*)pCur;
	//pCur += iSize;

	//iSize = MAX_KEY_POINT_COUNT * sizeof(Sift_Keypoint);
	//iSize = ALIGN_SIZE_128(iSize);
	//oSift.m_pKey_Point = (Sift_Keypoint*)pCur;
	//pCur += iSize;
	//oSift.m_pBuffer_Start = pCur;

	//oSift.m_pBuffer_Start = pBuffer;
	unsigned int iSift_Size = iGet_Sift_Detect_Size(iWidth, iHeight, o_min);
	oSift.m_pBuffer_Start = (unsigned char*)pMalloc(iSift_Size);
	*poSift = oSift;
}

void Arrange_Mem(Sift* poSift, int w, int h, unsigned char* pCur)
{//对本轮进行 octave, DoG, grad进行内存安排
	Sift oSift = *poSift;
	int i, iSize, iSize_1;
	//unsigned char* pOrg = pCur;

	iSize = ALIGN_SIZE_128(w * h * sizeof(float));
	iSize_1 = ALIGN_SIZE_128(iSize + w * sizeof(float));

	//高斯图与DoG梅花间竹,尝试以下内存安排
	//|  /   | G_0 | DoG_0 | G_1 | DoG_1 | G_2
	//|DoG_2 | G_3 | DoG_3 | G_4 | DoG_4 | G_5
	if (oSift.o_cur == oSift.o_min)
	{
		for (i = 0; i < 3; i++)
		{
			oSift.m_pDoG[(i + 2) % 3] = (float*)pCur;
			pCur += iSize_1;
			oSift.m_pGauss[i] = (float*)pCur;
			pCur += iSize;
		}
	}
	else
	{//原来的[0]图做了一次的下采样,放在了[1],所以，循环左移一次
		float* pTemp[2] = { oSift.m_pGauss[0],oSift.m_pDoG[0] };
		for (i = 1; i < 3; i++)
		{
			oSift.m_pGauss[i - 1] = oSift.m_pGauss[i];
			oSift.m_pDoG[i - 1] = oSift.m_pDoG[i];
		}
		oSift.m_pGauss[2] = pTemp[0];
		oSift.m_pDoG[2] = pTemp[1];
	}
	*poSift = oSift;
}
#define Get_w_h(iWidth, iHeight,  w,  h, o_cur) \
{ \
	if (o_cur < 0) \
	{	\
		w = iWidth << (-o_cur);		\
		h = iHeight << (-o_cur);	\
	}else \
	{ \
		w = iWidth >> o_cur; \
		h = iHeight >> o_cur; \
	} \
}

static void Resize_Image_x2(float* pSource, int iWidth, int iHeight, float* pDest)
{//图像按边长扩大一倍
	//先行扩展
	int y, x, iWidth_x2 = iWidth << 1, iHeight_x2 = iHeight << 1,
		iHeight_x2_Minus_1 = iHeight_x2 - 1,
		iWidth_Minus_1 = iWidth - 1;
	float* pSource_Cur, * pDest_Cur;
	for (pSource_Cur = pSource, y = 0; y < iHeight; y++/*pSource_Cur+=iWidth*/)
	{
		pDest_Cur = &pDest[(y << 1) * iWidth_x2];
		for (x = 0; x < iWidth_Minus_1; x++, pSource_Cur++, pDest_Cur += 2)
		{
			pDest_Cur[0] = pSource_Cur[0];
			pDest_Cur[1] = (pSource_Cur[0] + pSource_Cur[1]) * 0.5f;
		}
		pDest_Cur[0] = pDest_Cur[1] = *pSource_Cur++;
	}

	//再对pDest再做一次的插值
	pDest_Cur = &pDest[iWidth_x2];
	for (y = 1; y < iHeight_x2_Minus_1; y += 2)
	{
		for (x = 0; x < iWidth_x2; x++, pDest_Cur++)
			pDest_Cur[0] = (pDest_Cur[-iWidth_x2] + pDest_Cur[iWidth_x2]) * 0.5f;
		pDest_Cur += iWidth_x2;
	}
	for (x = 0; x < iWidth_x2; x++, pDest_Cur++)
		*pDest_Cur = pDest_Cur[-iWidth_x2];
	return;
}
static void Resize_Image_Half(float* pSource, int iWidth, int iHeight, float* pDest)
{//边长减半
	int iWidth_Half = iWidth >> 1,
		iHeight_Half = iHeight >> 1;
	int y;
	float* pSource_Cur, * pDest_Cur = pDest, * pDest_End;
	for (y = 0; y < iHeight_Half; y++)
	{
		pSource_Cur = &pSource[y * iWidth << 1];
		pDest_End = pDest_Cur + iWidth_Half;
		for (; pDest_Cur < pDest_End; pDest_Cur++, pSource_Cur += 2)
			*pDest_Cur = *pSource_Cur;
	}
	return;
}

static void Sift_Dectact(Sift* poSift, int iCur_Layer)
{//这个和原来的Sift_Detect不同，做三层DoG
	Sift oSift = *poSift;
	int i, iSize;

	//生成一个sigma
	float fSigma_a, fSigma_b;
	//从以下代码看出高斯核的大小与图像大小无关
	if (oSift.o_cur == oSift.o_min)
	{//相当于Process_First_Octave
		fSigma_a = oSift.sigma0 * (float)pow(oSift.sigmak, oSift.s_min);
		fSigma_b = oSift.sigman * (float)pow(2.0f, -oSift.o_min);
	}
	else
	{//相当于Process_Next_Octave
		int s_best = min(oSift.s_min + oSift.S, oSift.s_max);
		fSigma_a = (float)(oSift.sigma0 * pow((float)oSift.sigmak, (float)oSift.s_min));
		fSigma_b = (float)(oSift.sigma0 * pow((float)oSift.sigmak, (float)(s_best - oSift.S)));
	}
	float sd;
	int r;
	float* pFilter;

	//先搞第0层Octave，重用这一层对应的DoG
	oSift.m_pTemp = oSift.m_pDoG[0];
	if (fSigma_a > fSigma_b)
	{
		sd = (float)sqrt(fSigma_a * fSigma_a - fSigma_b * fSigma_b);
		//可见此处是r随sigma变化而变化
		r = (int)max((float)ceil(4.0 * sd), 1.f);
		Gen_Gauss_Filter(r, sd, &pFilter);
		Gauss_Filter_AVX512(oSift.m_pGauss[0], oSift.octave_width, oSift.octave_height, oSift.m_pGauss[0], oSift.m_pTemp, r, pFilter);
		Free(pFilter);
	}
	union {
		struct {
			float* pt, * pOctave_0, * pOctave_1;
		};
		__m256i Addr_3;
	};
	float* pEnd;
	iSize = oSift.octave_width * oSift.octave_height;

	for (i = 0; i < 2; i++)
	{//连做两层的Octave和DoG
		pOctave_0 = oSift.m_pGauss[i];
		pOctave_1 = oSift.m_pGauss[1 + i];
		oSift.m_pTemp = oSift.m_pDoG[1 + i];
		sd = oSift.dsigma0 * (float)pow(oSift.sigmak, i);

		r = (int)max((float)ceil(4.0 * sd), 1.f);
		Gen_Gauss_Filter(r, sd, &pFilter);
		//从下面一层生成上面一层
		Gauss_Filter_AVX512(pOctave_0, oSift.octave_width, oSift.octave_height, pOctave_1, oSift.m_pTemp, r, pFilter);
		Free(pFilter);

		//搞第0层DoG
		pt = oSift.m_pDoG[i];
		pEnd = pt + iSize;
		for (; pt < pEnd; Addr_3 = _mm256_add_epi64(Addr_3, _mm256_set1_epi64x(16 * 4))   /*pt += 16, pOctave_0 += 16, pOctave_1 += 16*/)
			*(__m512*)pt = _mm512_sub_ps(*(__m512*)pOctave_1, *(__m512*)pOctave_0);
	}
}
static int Resize_Keypoint_Feature_Size(Sift* poSift, Sift_Feature_Level* poFeature_Level,int bKeypoint)
{//增加分配内存, bKeypoint =1 给Keypoint 加内存	=0 给Feature 加
	int iResult = 0, iOffset;
	if (bKeypoint)
	{
		Sift_Keypoint* pOrg = poFeature_Level->m_pKeypoint;
		iOffset = (int)(pOrg - poSift->m_pKey_Point);
		poSift->m_iMax_Keypoint_Count *= 2;
		//Disp_Mem();
		iResult = bExpand_Ex((void**)&poSift->m_pKey_Point, poSift->m_iMax_Keypoint_Count * sizeof(Sift_Keypoint));
		if (iResult)
			poFeature_Level->m_pKeypoint = poSift->m_pKey_Point + iOffset;
		else
			printf("Fail to expamd memroy in Resize_Keypoint_Feature_Size\n");
	}else
	{
		Sift_Feature* pOrg = poFeature_Level->m_pFeature;
		iOffset = (int)( pOrg - poSift->m_pFeature);
		poSift->m_iMax_Feature_Count *= 2;
		iResult = bExpand_Ex((void**)&poSift->m_pFeature, poSift->m_iMax_Feature_Count * sizeof(Sift_Feature));
		if (iResult)
			poFeature_Level->m_pFeature = poSift->m_pFeature + iOffset;
	}
	return iResult;
}
static int Get_Keypoint(Sift* poSift, int iDoG_Mid, int s, Sift_Feature_Level* poFeature_Level)
{//注意： iDoG_Mid指的是本次处理DoG中间层在oSift.m_pDoG中逻辑位置
//			s:对应sift源代码中的s
	Sift oSift = *poSift;
	int iSize = oSift.octave_width * oSift.octave_height;
	int const    xo = 1;					/* x-stride */
	int const    yo = oSift.octave_width;	/* y-stride */
	iSize = ALIGN_SIZE_128(iSize) * 2;
	//此处可能要改，不能用原来的方法
	int iDoG_Mid_Index = iDoG_Mid % 3;
	int so_Next;	// = iDoG_Mid_Index == 2 ? -2 * iSize : iSize;
	int so_Pre;		// = iDoG_Mid_Index == 0 ? iSize * 2 : -iSize;
	float tp = oSift.peak_thresh * 0.8f;
	int x_End_Align_16 = 1 + (((oSift.octave_width - 2) >> 4) << 4),
		iMask = (1 << ((oSift.octave_width - 2) & 0xF)) - 1;
	int y, x;
	int i, j, k, ii, jj;
	float* pt;

	so_Next = (int)(oSift.m_pDoG[(iDoG_Mid + 1) % 3] - oSift.m_pDoG[iDoG_Mid_Index]);
	so_Pre = (int)(oSift.m_pDoG[(iDoG_Mid - 1) % 3] - oSift.m_pDoG[iDoG_Mid_Index]);

	//以DoG[1]层为中心层，通过26邻域确定最初的Key_Point
	pt = &oSift.m_pDoG[iDoG_Mid_Index][1 * oSift.octave_width + 1];
	Sift_Keypoint oKey_Point, * pPoint = poFeature_Level->m_pKeypoint;	//oSift.m_pKey_Point + oSift.nkeys;
	int nkeys = 0;
	int iRemain = oSift.m_iMax_Keypoint_Count - (int)(pPoint - oSift.m_pKey_Point);	//还剩下读诵好个
	int Offset_26[] = { -xo, xo,
						-yo - xo, -yo, -yo + xo,
						yo - xo,	yo, yo + xo,

						-xo + so_Pre, so_Pre, xo + so_Pre,
						-yo - xo + so_Pre, -yo + so_Pre, -yo + xo + so_Pre,
						yo - xo + so_Pre,	yo + so_Pre, yo + xo + so_Pre,

						-xo + so_Next, so_Next, xo + so_Next,
						-yo - xo + so_Next, -yo + so_Next, -yo + xo + so_Next,
						yo - xo + so_Next,	yo + so_Next, yo + xo + so_Next
	};

	__m512 v_16;

#define CHECK_NEIGHBORS_1(CMP,SGN)   (    \
_mm512_cmp_ps_mask(v_16, _mm512_loadu_ps(pt + Offset_26[0]),CMP) && \
_mm512_cmp_ps_mask(v_16, _mm512_loadu_ps(pt + Offset_26[1]),CMP) && \
_mm512_cmp_ps_mask(v_16, _mm512_loadu_ps(pt + Offset_26[2]),CMP) && \
_mm512_cmp_ps_mask(v_16, _mm512_loadu_ps(pt + Offset_26[3]),CMP) && \
_mm512_cmp_ps_mask(v_16, _mm512_loadu_ps(pt + Offset_26[4]),CMP) && \
_mm512_cmp_ps_mask(v_16, _mm512_loadu_ps(pt + Offset_26[5]),CMP) && \
_mm512_cmp_ps_mask(v_16, _mm512_loadu_ps(pt + Offset_26[6]),CMP) && \
_mm512_cmp_ps_mask(v_16, _mm512_loadu_ps(pt + Offset_26[7]),CMP) && \
\
_mm512_cmp_ps_mask(v_16, _mm512_loadu_ps(pt + Offset_26[8]),CMP) && \
_mm512_cmp_ps_mask(v_16, _mm512_loadu_ps(pt + Offset_26[9]),CMP) && \
_mm512_cmp_ps_mask(v_16, _mm512_loadu_ps(pt + Offset_26[10]),CMP) && \
_mm512_cmp_ps_mask(v_16, _mm512_loadu_ps(pt + Offset_26[11]),CMP) && \
_mm512_cmp_ps_mask(v_16, _mm512_loadu_ps(pt + Offset_26[12]),CMP) && \
_mm512_cmp_ps_mask(v_16, _mm512_loadu_ps(pt + Offset_26[13]),CMP) && \
_mm512_cmp_ps_mask(v_16, _mm512_loadu_ps(pt + Offset_26[14]),CMP) && \
_mm512_cmp_ps_mask(v_16, _mm512_loadu_ps(pt + Offset_26[15]),CMP) && \
_mm512_cmp_ps_mask(v_16, _mm512_loadu_ps(pt + Offset_26[16]),CMP)  &&\
\
_mm512_cmp_ps_mask(v_16, _mm512_loadu_ps(pt + Offset_26[17]),CMP) && \
_mm512_cmp_ps_mask(v_16, _mm512_loadu_ps(pt + Offset_26[18]),CMP) && \
_mm512_cmp_ps_mask(v_16, _mm512_loadu_ps(pt + Offset_26[19]),CMP) && \
_mm512_cmp_ps_mask(v_16, _mm512_loadu_ps(pt + Offset_26[20]),CMP) && \
_mm512_cmp_ps_mask(v_16, _mm512_loadu_ps(pt + Offset_26[21]),CMP) && \
_mm512_cmp_ps_mask(v_16, _mm512_loadu_ps(pt + Offset_26[22]),CMP) && \
_mm512_cmp_ps_mask(v_16, _mm512_loadu_ps(pt + Offset_26[23]),CMP) && \
_mm512_cmp_ps_mask(v_16, _mm512_loadu_ps(pt + Offset_26[24]),CMP) && \
_mm512_cmp_ps_mask(v_16, _mm512_loadu_ps(pt + Offset_26[25]),CMP) && \
_mm512_cmp_ps_mask(v_16, _mm512_set1_ps(SGN tp),CMP)  \
		) \

#define CHECK_NEIGHBORS(CMP,SGN)   (    \
_mm512_cmp_ps_mask(v_16, _mm512_set1_ps(SGN tp),CMP) & \
_mm512_cmp_ps_mask(v_16, _mm512_loadu_ps(pt + Offset_26[0]),CMP) & \
_mm512_cmp_ps_mask(v_16, _mm512_loadu_ps(pt + Offset_26[1]),CMP) & \
_mm512_cmp_ps_mask(v_16, _mm512_loadu_ps(pt + Offset_26[2]),CMP) & \
_mm512_cmp_ps_mask(v_16, _mm512_loadu_ps(pt + Offset_26[3]),CMP) & \
_mm512_cmp_ps_mask(v_16, _mm512_loadu_ps(pt + Offset_26[4]),CMP) & \
_mm512_cmp_ps_mask(v_16, _mm512_loadu_ps(pt + Offset_26[5]),CMP) & \
_mm512_cmp_ps_mask(v_16, _mm512_loadu_ps(pt + Offset_26[6]),CMP) & \
_mm512_cmp_ps_mask(v_16, _mm512_loadu_ps(pt + Offset_26[7]),CMP) & \
\
_mm512_cmp_ps_mask(v_16, _mm512_loadu_ps(pt + Offset_26[8]),CMP) & \
_mm512_cmp_ps_mask(v_16, _mm512_loadu_ps(pt + Offset_26[9]),CMP) & \
_mm512_cmp_ps_mask(v_16, _mm512_loadu_ps(pt + Offset_26[10]),CMP) & \
_mm512_cmp_ps_mask(v_16, _mm512_loadu_ps(pt + Offset_26[11]),CMP) & \
_mm512_cmp_ps_mask(v_16, _mm512_loadu_ps(pt + Offset_26[12]),CMP) & \
_mm512_cmp_ps_mask(v_16, _mm512_loadu_ps(pt + Offset_26[13]),CMP) & \
_mm512_cmp_ps_mask(v_16, _mm512_loadu_ps(pt + Offset_26[14]),CMP) & \
_mm512_cmp_ps_mask(v_16, _mm512_loadu_ps(pt + Offset_26[15]),CMP) & \
_mm512_cmp_ps_mask(v_16, _mm512_loadu_ps(pt + Offset_26[16]),CMP)  &\
\
_mm512_cmp_ps_mask(v_16, _mm512_loadu_ps(pt + Offset_26[17]),CMP) & \
_mm512_cmp_ps_mask(v_16, _mm512_loadu_ps(pt + Offset_26[18]),CMP) & \
_mm512_cmp_ps_mask(v_16, _mm512_loadu_ps(pt + Offset_26[19]),CMP) & \
_mm512_cmp_ps_mask(v_16, _mm512_loadu_ps(pt + Offset_26[20]),CMP) & \
_mm512_cmp_ps_mask(v_16, _mm512_loadu_ps(pt + Offset_26[21]),CMP) & \
_mm512_cmp_ps_mask(v_16, _mm512_loadu_ps(pt + Offset_26[22]),CMP) & \
_mm512_cmp_ps_mask(v_16, _mm512_loadu_ps(pt + Offset_26[23]),CMP) & \
_mm512_cmp_ps_mask(v_16, _mm512_loadu_ps(pt + Offset_26[24]),CMP) & \
_mm512_cmp_ps_mask(v_16, _mm512_loadu_ps(pt + Offset_26[25]),CMP)  \
		) \

	int iResult;
	for (y = 1; y < oSift.octave_height - 1; y++)
	{
		pt = &oSift.m_pDoG[iDoG_Mid_Index][y * oSift.octave_width + 1];
		for (x = 1; x < x_End_Align_16; x += 16, pt += 16)
		{
			//v_16 = *(__m512*)pt;
			v_16 = _mm512_loadu_ps(pt);

			if (/*iResult =*/ CHECK_NEIGHBORS_1(_CMP_GT_OS, +) || CHECK_NEIGHBORS_1(_CMP_LT_OS, -))
			{
				iResult = CHECK_NEIGHBORS(_CMP_GT_OS, +) | CHECK_NEIGHBORS(_CMP_LT_OS, -);
				for (i = 0; iResult; i++)
				{
					if (iResult & 1)
					{
						oKey_Point.ix = x + i;
						oKey_Point.iy = y;
						oKey_Point.is = s;
						if (iRemain <= 0)
						{//需要扩大Size
							if (!(Resize_Keypoint_Feature_Size(&oSift, poFeature_Level, 1)))
								return 0;
							pPoint = poFeature_Level->m_pKeypoint;
							iRemain += oSift.m_iMax_Keypoint_Count >> 1;
						}
						iRemain--;
						pPoint[nkeys] = oKey_Point;
						nkeys++;
						
					}
					iResult >>= 1;
				}
			}
		}

		//搞余数
		//v_16 = *(__m512*)pt;
		v_16 = _mm512_loadu_ps(pt);
		iResult = (CHECK_NEIGHBORS(_CMP_GT_OS, +) | CHECK_NEIGHBORS(_CMP_LT_OS, -)) & iMask;
		for (i = 0; iResult; i++)
		{
			if (iResult & 1)
			{
				oKey_Point.ix = x + i;
				oKey_Point.iy = y;
				oKey_Point.is = s;
				if (iRemain <= 0)
				{//需要扩大Size
					if (!(Resize_Keypoint_Feature_Size(&oSift, poFeature_Level, 1)))
						return 0;
					pPoint = poFeature_Level->m_pKeypoint;
					iRemain += oSift.m_iMax_Keypoint_Count >> 1;
				}
				iRemain--;
				pPoint[nkeys] = oKey_Point;
				nkeys++;
			}
			iResult >>= 1;
		}
	}

	//对上一步的keypoints进一步筛选
	int dx, dy, iter;
	float Dx = 0, Dy = 0, Ds = 0, Dxx = 0, Dyy = 0, Dss = 0, Dxy = 0, Dxs = 0, Dys = 0;
	float A[3][3], b[3];	//正是求解Ax=b，可以直接调用解线性方程
	int iValid_Point_Count = 0;
	float* pDoG = oSift.m_pDoG[iDoG_Mid_Index];
	float xper = (float)pow(2.0, oSift.o_cur);
	iRemain = oSift.m_iMax_Keypoint_Count - (int)(pPoint - oSift.m_pKey_Point);	//还剩下读诵好个
	for (k = 0; k < nkeys; k++)
	{
		oKey_Point = pPoint[k];
		//先来个一阶导
		Dx = 0, Dy = 0, Ds = 0, Dxx = 0, Dyy = 0, Dss = 0, Dxy = 0, Dxs = 0, Dys = 0;
		dx = dy = 0;
		x = oKey_Point.ix;
		y = oKey_Point.iy;

		for (iter = 0; iter < 5; iter++)
		{//感觉迭代多余，也有可能为了不落在局部极小点或着噪声影响等等原因，再继续搜索一下
			//而且一点的x,y各自方向上最大跨度不可能大于1，所以在9邻域上游走
			x += dx;
			y += dy;
			pt = pDoG + xo * x + yo * y;

			//显然就是求一阶偏导，雅比克
			Dx = 0.5f * (pt[1] - pt[-1]);
			Dy = 0.5f * (pt[oSift.octave_width] - pt[-oSift.octave_width]);
			Ds = 0.5f * (pt[so_Next] - pt[so_Pre]);

			//求一个Hesse矩阵，由二阶偏导组成
			Dxx = (pt[1] + pt[-1] - 2.0f * (*pt));
			Dyy = (pt[oSift.octave_width] + pt[-oSift.octave_width] - 2.0f * (*pt));
			Dss = (pt[so_Next] + pt[so_Pre] - 2.0f * (*pt));

			Dxy = 0.25f * (pt[oSift.octave_width + 1] + pt[-oSift.octave_width - 1] - pt[-oSift.octave_width + 1] - pt[oSift.octave_width - 1]);	//  at(+1, +1, 0) + at(-1, -1, 0) - at(-1, +1, 0) - at(+1, -1, 0));
			Dxs = 0.25f * (pt[so_Next + 1] + pt[so_Pre - 1] - pt[so_Pre + 1] - pt[so_Next - 1]);	//  (at(+1, 0, +1) + at(-1, 0, -1) - at(-1, 0, +1) - at(+1, 0, -1));
			Dys = 0.25f * (pt[so_Next + oSift.octave_width] + pt[so_Pre - oSift.octave_width] - pt[so_Pre + oSift.octave_width] - pt[so_Next - oSift.octave_width]);	//at(0, +1, +1) + at(0, -1, -1) - at(0, -1, +1) - at(0, +1, -1));

			//搞个对称阵hesse
			A[0][0] = Dxx;
			A[1][1] = Dyy;
			A[2][2] = Dss;
			A[0][1] = A[1][0] = Dxy;
			A[0][2] = A[2][0] = Dxs;
			A[1][2] = A[2][1] = Dys;

			b[0] = -Dx;
			b[1] = -Dy;
			b[2] = -Ds;

			//此处为泰勒修正。以上找到的关键点并非真正的极值点，在插值的条件下，极值点
			//应该是子像素，也就是落在以上关键点邻域上。故此要找在关键点周围做一个泰勒展开
			//得到一个近似 f(x)= f(x0) + J*(delta x) + 1/2(delta x') * H * (delta x)
			//对 f(x) 求极值 min f(x) 就是再来一次求偏导 J + H*(delta x) ，使其=0, 求 delta x
			//所以有 H*x= -J ,解出 x就是所求的 delta x
			for (j = 0; j < 3; j++)
			{//j为列
				float maxa = 0;
				float maxabsa = 0;
				int    maxi = -1;
				float tmp;

				/* look for the maximally stable pivot */
				for (i = j; i < 3; i++)
				{//显然找第j列最大绝对值，列主元
					float a = A[i][j];
					float absa = abs(a);
					if (absa > maxabsa)
					{
						maxa = a;
						maxabsa = absa;
						maxi = i;
					}
				}

				//A为二阶导矩阵，如果二阶导矩阵存在一列<eps则退出？表示这是这个方向上存在极值点?
				if (maxabsa < 1e-10f)
				{
					b[0] = 0;
					b[1] = 0;
					b[2] = 0;
					break;
				}
				i = maxi;

				//这里越看越象列主元法解线性方程 Ax=b
				/* swap j-th row with i-th row and normalize j-th row */
				for (jj = j; jj < 3; ++jj)
				{//将最大元素所在行与第j行交换，并且规格化
					tmp = A[i][jj]; A[i][jj] = A[j][jj]; A[j][jj] = tmp;
					A[j][jj] /= maxa;
				}
				tmp = b[j]; b[j] = b[i]; b[i] = tmp;
				b[j] /= maxa;

				for (ii = j + 1; ii < 3; ++ii)
				{
					float x = A[ii][j];
					for (jj = j; jj < 3; ++jj)
						A[ii][jj] -= x * A[j][jj];
					b[ii] -= x * b[j];
				}
			}

			// backward substitution 回代，瞅着眼熟
			for (i = 2; i > 0; --i)
			{
				float x = b[i];
				for (ii = i - 1; ii >= 0; --ii)
					b[ii] -= x * A[ii][i];
			}

			dx = ((b[0] > 0.6 && x < oSift.octave_width - 2) ? 1 : 0)
				+ ((b[0] < -0.6 && x > 1) ? -1 : 0);

			dy = ((b[1] > 0.6 && y < oSift.octave_height - 2) ? 1 : 0)
				+ ((b[1] < -0.6 && y > 1) ? -1 : 0);

			if (dx == 0 && dy == 0)
				break;
		}

		/* check threshold and other conditions */
		{
			float val = *pt + 0.5f * (Dx * b[0] + Dy * b[1] + Ds * b[2]);
			//这句的意义待分析
			float score = (Dxx + Dyy) * (Dxx + Dyy) / (Dxx * Dyy - Dxy * Dxy);
			float xn = x + b[0];
			float yn = y + b[1];
			float sn = s + b[2];
			// |H|=Dxx * Dyy - Dxy * Dxy ,

			int good = abs(val) > oSift.peak_thresh &&
				score < (oSift.edge_thresh + 1) * (oSift.edge_thresh + 1) / oSift.edge_thresh &&
				score >= 0 &&
				abs(b[0]) < 1.5 &&
				abs(b[1]) < 1.5 &&
				abs(b[2]) < 1.5 &&
				xn >= 0 &&
				xn <= oSift.octave_width - 1 &&
				yn >= 0 &&
				yn <= oSift.octave_height - 1 &&
				sn >= oSift.s_min &&
				sn <= oSift.s_max;

			if (good)
			{
				oKey_Point.ix = x;
				oKey_Point.iy = y;
				oKey_Point.is = s;
				oKey_Point.o = oSift.o_cur;
				oKey_Point.s = sn;
				oKey_Point.x = xn * xper;
				oKey_Point.y = yn * xper;
				oKey_Point.sigma = oSift.sigma0 * (float)pow(2.0, sn / oSift.S) * xper;
				//if (iRemain <= 0)
				//{//需要扩大Size
				//	if (!(Resize_Keypoint_Feature_Size(&oSift, poFeature_Level, 1)))
				//		return 0;
				//	pPoint = poFeature_Level->m_pKeypoint;
				//	iRemain += oSift.m_iMax_Keypoint_Count >> 1;
				//}
				//iRemain--;
				pPoint[iValid_Point_Count++] = oKey_Point;
			}
		} /* done checking */
	}

	poFeature_Level->m_iKeypoint_Count = iValid_Point_Count;
	*poSift = oSift;

	return 1;
#undef CHECK_NEIGHBORS
}

void _Update_Gradiant_AVX512(float* src, float* grad, int w, int h)
{//个人认为此算法可以更简洁，从sift的特征点看，角点的意义最大，边缘点和孤立点都丢弃，故此落在
	//图像上下四边上的点基本上不会是角点。如果后面有充分的数据验证这个猜想，则边上的数值可以不算
	float* pSource, * pSource_1, * pEnd_1, * pEnd/*,*pEnd_Align_16*/;
	float* pGrad, * pGrad_1;
	__m512 gx, gy, Value;
	__m512i vIndex_Grad;
	int i, iResult, iValue, h_Minus_1 = h - 1, w_x2 = w * 2;
	for (i = 0; i < 16; i++)
		((int*)&vIndex_Grad)[i] = i * 8;
	//vIndex_Grad.m512i_i32[i] = i * 8;

//先不管上下左右角4个点，有必要省略掉

//第一点，没什么好搞了
/*gx.m512_f32[0] = src[1] - src[0];
gy.m512_f32[0] = src[w] - src[0];
grad[0] = (float)sqrt(gx.m512_f32[0] * gx.m512_f32[0] + gy.m512_f32[0] * gy.m512_f32[0]);
grad[1] = (float)atan(gy.m512_f32[0] / gx.m512_f32[0]);*/
	float fgx = src[1] - src[0],
		fgy = src[w] - src[0];
	grad[0] = (float)sqrt(fgx * fgx + fgy * fgy);
	grad[1] = (float)atan(fgy / fgx);
	if (grad[1] < 0)
		grad[1] +=(float)( PI * 2.f);

	//搞第一行
	pGrad = grad + 2;
	pSource = src + 1;
	pEnd = src + w - 1;	//头尾两点不搞
	//pEnd_Align_16 = pSource + (unsigned int)((pEnd - pSource) & 0xFFFFFFF0);
	iValue = (pEnd - pSource) & 0xF;
	while (pSource < pEnd)
	{//放心推进，有余数无所谓，留有余地即可
		//gx = 0.5f * (src[+xo] - src[-xo]);
		gx = _mm512_mul_ps(_mm512_set1_ps(0.5f), _mm512_sub_ps(_mm512_loadu_ps(pSource + 1), *(__m512*)(pSource - 1)));
		//gy = src[+yo] - src[0];
		gy = _mm512_sub_ps(_mm512_loadu_ps(pSource + w), _mm512_loadu_ps(pSource));

		//sqrt(gx * gx + gy * gy)
		Value = _mm512_sqrt_ps(_mm512_add_ps(_mm512_mul_ps(gx, gx), _mm512_mul_ps(gy, gy)));
		_mm512_i32scatter_ps(pGrad, vIndex_Grad, Value, 1);

		//此处感觉可以简化 arctg的判断，值域为(-PI/2,PI/2);
		Value = _mm512_atan_ps(_mm512_div_ps(gy, gx));
		iResult = _mm512_cmp_ps_mask(Value, _mm512_set1_ps(0), _CMP_LT_OS);
		Value = _mm512_mask_add_ps(Value, iResult, Value, _mm512_set1_ps((float)(PI * 2)));
		_mm512_i32scatter_ps(pGrad + 1, vIndex_Grad, Value, 1);
		pSource += 16;
		pGrad += 32;
	}

	//第一行最后一点
	//gx = src[0] - src[-xo];
	//gy = src[+yo] - src[0];
	pSource = &src[w - 1];
	pGrad = &grad[w_x2 - 2];
	/*gx.m512_f32[0] = pSource[0] - pSource[-1];
	gy.m512_f32[0] = pSource[w] - pSource[0];
	pGrad[0] = (float)sqrt(gx.m512_f32[0] * gx.m512_f32[0] + gy.m512_f32[0] * gy.m512_f32[0]);
	pGrad[1] = (float)atan(gy.m512_f32[0] / gx.m512_f32[0]);*/
	fgx = pSource[0] - pSource[-1];
	fgy = pSource[w] - pSource[0];
	pGrad[1] = (float)sqrt(fgx * fgx + fgy * fgy);
	pGrad[1] = (float)atan(fgy / fgx);

	if (pGrad[1] < 0)
		pGrad[1] += (float)(PI * 2.f);

	//搞中间行
	pSource_1 = src + w + 1;
	pEnd_1 = src + h_Minus_1 * w;
	pGrad_1 = grad + w * 2 + 2;
	pEnd = pSource_1 + w - 2;
	//以下代码还有Bug, 以下循环中的最后一行还是会出现踩过界算错问题
	for (; pSource_1 < pEnd_1; pSource_1 += w, pGrad_1 += w_x2, pEnd += w)
	{
		pSource = pSource_1;
		pGrad = pGrad_1;

		//中间行第一点
		//gx = src[+xo] - src[0];
		//gy = 0.5f * (src[+yo] - src[-yo]);
		/*gx.m512_f32[0] = pSource[0] - pSource[-1];
		gy.m512_f32[0] = 0.5f * (pSource[w - 1] - pSource[-1 - w]);
		pGrad[-2] = (float)sqrt(gx.m512_f32[0] * gx.m512_f32[0] + gy.m512_f32[0] * gy.m512_f32[0]);
		pGrad[-1] = (float)atan(gy.m512_f32[0] / gx.m512_f32[0]);*/
		fgx = pSource[0] - pSource[-1];
		fgy = 0.5f * (pSource[w - 1] - pSource[-1 - w]);
		pGrad[-2] = (float)sqrt(fgx * fgx + fgy * fgy);
		pGrad[-1] = (float)atan(fgy / fgx);
		if (pGrad[-1] < 0)
			pGrad[-1] +=(float)( PI * 2.f);

		while (pSource < pEnd)
		{
			//gx = 0.5f * (src[+xo] - src[-xo]);	//根本没有按3x3模板搞，偷工减料做了一个一维的模板
			gx = _mm512_mul_ps(_mm512_set1_ps(0.5f), _mm512_sub_ps(_mm512_loadu_ps(pSource + 1), _mm512_loadu_ps(pSource - 1)));
			//gy = 0.5f * (src[+yo] - src[-yo]);	//同上
			gy = _mm512_mul_ps(_mm512_set1_ps(0.5f), _mm512_sub_ps(_mm512_loadu_ps(pSource + w), _mm512_loadu_ps(pSource - w)));

			//sqrt(gx * gx + gy * gy)
			Value = _mm512_sqrt_ps(_mm512_add_ps(_mm512_mul_ps(gx, gx), _mm512_mul_ps(gy, gy)));
			_mm512_i32scatter_ps(pGrad, vIndex_Grad, Value, 1);

			//此处感觉可以简化 arctg的判断，值域为(-PI/2,PI/2);
			Value = _mm512_atan_ps(_mm512_div_ps(gy, gx));
			iResult = _mm512_cmp_ps_mask(Value, _mm512_set1_ps(0), _CMP_LE_OS);
			Value = _mm512_mask_add_ps(Value, iResult, Value, _mm512_set1_ps((float)(PI * 2)));
			_mm512_i32scatter_ps(pGrad + 1, vIndex_Grad, Value, 1);

			pSource += 16;
			pGrad += 32;
		}

		//中间行最后一点
		pSource = &pSource_1[w - 2];
		pGrad = &pGrad_1[w_x2 - 4];
		//gx = src[0] - src[-xo];
		//gy = 0.5f * (src[+yo] - src[-yo]);
		/*gx.m512_f32[0] = pSource[0] - pSource[-1];
		gy.m512_f32[0] = 0.5f * (pSource[w] - pSource[-w]);
		pGrad[0] = (float)sqrt(gx.m512_f32[0] * gx.m512_f32[0] + gy.m512_f32[0] * gy.m512_f32[0]);
		pGrad[1] = (float)atan(gy.m512_f32[0] / gx.m512_f32[0]);*/
		fgx = pSource[0] - pSource[-1];
		fgy = 0.5f * (pSource[w] - pSource[-w]);
		pGrad[0] = (float)sqrt(fgx * fgx + fgy * fgy);
		pGrad[1] = (float)atan(fgy / fgx);

		if (pGrad[1] < 0)
			pGrad[1] +=(float)( PI * 2.f);
	}

	//最后一行
	pGrad = pGrad_1;
	pSource = pSource_1;
	//pEnd = pSource_1 + w - 2;	//头尾两点不搞
	iValue = (((w - 2) >> 4) << 4);
	pEnd = pSource + iValue;
	iValue = (1 << (w - 2 - iValue)) - 1;

	//最后一行第一点
	//gx = src[+xo] - src[0];
	//gy = src[0] - src[-yo];
	/*gx.m512_f32[0] = pSource[0] - pSource[-1];
	gy.m512_f32[0] = pSource[-1] - pSource[-1 - w];
	pGrad[-2] = (float)sqrt(gx.m512_f32[0] * gx.m512_f32[0] + gy.m512_f32[0] * gy.m512_f32[0]);
	pGrad[-1] = (float)atan(gy.m512_f32[0] / gx.m512_f32[0]);*/
	fgx = pSource[0] - pSource[-1];
	fgy = pSource[-1] - pSource[-1 - w];
	pGrad[-2] = (float)sqrt(fgx * fgx + fgy * fgy);
	pGrad[-1] = (float)atan(fgy / fgx);
	if (pGrad[-1] < 0)
		pGrad[-1] += (float)(PI * 2.f);

	while (pSource < pEnd)
	{
		//gx = 0.5f * (src[+xo] - src[-xo]);
		gx = _mm512_mul_ps(_mm512_set1_ps(0.5f), _mm512_sub_ps(_mm512_loadu_ps(pSource + 1), _mm512_loadu_ps(pSource - 1)));
		//gy = src[0] - src[-yo];
		gy = _mm512_sub_ps(_mm512_loadu_ps(pSource), _mm512_loadu_ps(pSource - w));

		//sqrt(gx * gx + gy * gy)
		Value = _mm512_sqrt_ps(_mm512_add_ps(_mm512_mul_ps(gx, gx), _mm512_mul_ps(gy, gy)));
		_mm512_i32scatter_ps(pGrad, vIndex_Grad, Value, 1);

		//此处感觉可以简化 arctg的判断，值域为(-PI/2,PI/2);
		Value = _mm512_atan_ps(_mm512_div_ps(gy, gx));
		iResult = _mm512_cmp_ps_mask(Value, _mm512_set1_ps(0), _CMP_LT_OS);
		Value = _mm512_mask_add_ps(Value, iResult, Value, _mm512_set1_ps((float)(PI * 2)));
		_mm512_i32scatter_ps(pGrad + 1, vIndex_Grad, Value, 1);

		pSource += 16;
		pGrad += 32;
	}

	//最后一行余数
	gx = _mm512_mul_ps(_mm512_set1_ps(0.5f), _mm512_sub_ps(_mm512_loadu_ps(pSource + 1), _mm512_loadu_ps(pSource - 1)));
	gy = _mm512_sub_ps(_mm512_loadu_ps(pSource), _mm512_loadu_ps(pSource - w));

	Value = _mm512_sqrt_ps(_mm512_add_ps(_mm512_mul_ps(gx, gx), _mm512_mul_ps(gy, gy)));
	_mm512_mask_i32scatter_ps(pGrad, iValue, vIndex_Grad, Value, 1);

	//此处感觉可以简化 arctg的判断，值域为(-PI/2,PI/2);
	Value = _mm512_atan_ps(_mm512_div_ps(gy, gx));
	iResult = _mm512_cmp_ps_mask(Value, _mm512_set1_ps(0), _CMP_LT_OS);
	Value = _mm512_mask_add_ps(Value, iResult, Value, _mm512_set1_ps((float)(PI * 2)));
	_mm512_mask_i32scatter_ps(pGrad + 1, iValue, vIndex_Grad, Value, 1);

	//最后一行最后一点
	pSource = &pSource_1[w - 2];
	pGrad = &pGrad_1[w_x2 - 4];
	/*gx.m512_f32[0] = pSource[0] - pSource[-1];
	gy.m512_f32[0] = pSource[0] - pSource[-w];
	pGrad[0] = (float)sqrt(gx.m512_f32[0] * gx.m512_f32[0] + gy.m512_f32[0] * gy.m512_f32[0]);
	pGrad[1] = (float)atan(gy.m512_f32[0] / gx.m512_f32[0]);*/
	fgx = pSource[0] - pSource[-1];
	fgy = pSource[0] - pSource[-w];
	pGrad[0] = (float)sqrt(fgx * fgx + fgy * fgy);
	pGrad[1] = (float)atan(fgy / fgx);

	if (pGrad[-1] < 0)
		pGrad[-1] += (float)(PI * 2.f);

	////验算
	//FILE* pFile = fopen("c:\\tmp\\1.bin", "rb");
	//for (i = 0; i < w * h; i++)
	//{
	//	float Data[2];
	//	if (!fread(Data, 1, 2 * sizeof(float), pFile))
	//	{
	//		printf("error");
	//		break;
	//	}
	//	if (Data[0] != grad[i * 2] || Data[1] != grad[i * 2 + 1])
	//		printf("y:%d x:%d %f %f\n",i/w,i%w, Data[1], grad[i*2+1]);
	//}
	return;
}

int iCal_Key_Point_Oriatation(Sift* poSift, float Angle[4], Sift_Keypoint oKey_Point)
{//感觉此处是关键了
	float const winf = 1.5;
	float       xper = (float)pow(2.0, poSift->o_cur);

	int          w = poSift->octave_width;	//poSift->m_iWidth_x2;
	int          h = poSift->octave_height;	// poSift->m_iHeight_x2;
	int const    xo = 2;         /* x-stride */
	int const    yo = 2 * w;     /* y-stride */
	//int const    so = 2 * w * h; /* s-stride */
	float       x = oKey_Point.x / xper;
	float       y = oKey_Point.y / xper;
	float       sigma = oKey_Point.sigma / xper;

	int          xi = (int)(x + 0.5);
	int          yi = (int)(y + 0.5);
	int          si = oKey_Point.is;

	float const sigmaw = winf * sigma;
	int          W = (int)max((float)floor(3.0 * sigmaw), 1.f);
	int          nangles = 0;

	//enum { nbins = 36 };
	const int nbins = 36;
	float hist[nbins] = { 0 }, maxh;
	float const* pt;
	int xs, ys, iter, i;

	float dx;
	float dy;
	float r2;
	float wgt, mod, ang, fbin;
	/* skip if the keypoint octave is not current */
	if (oKey_Point.o != poSift->o_cur)
		return 0;

	//万一越界
	if (xi < 0 || xi > w - 1 || yi < 0 || yi > h - 1 ||
		si < poSift->s_min + 1 || si > poSift->s_max - 2)
		return 0;	//这个很奇怪，为什么会出现越界情况，感觉是伪判断

	/* clear histogram */
	//memset(hist, 0, sizeof(float) * nbins);

	/* compute orientation histogram */
	pt = poSift->m_pGrad + xo * xi + yo * yi;

#undef  at
#define at(dx,dy) (*(pt + xo * (dx) + yo * (dy)))

	for (ys = max(-W, -yi); ys <= min(+W, h - 1 - yi); ++ys)
	{
		for (xs = max(-W, -xi); xs <= min(+W, w - 1 - xi); ++xs)
		{
			dx = (float)(xi + xs) - x;
			dy = (float)(yi + ys) - y;
			r2 = dx * dx + dy * dy;
			//wgt, mod, ang, fbin;

			/* limit to a circular window */
			if (r2 >= W * W + 0.6)
				continue;

			//wgt = fast_expn(r2 / (2 * sigmaw * sigmaw));
			wgt = (float)exp(-r2 / (2 * sigmaw * sigmaw));
			//预先已经算出来，在Update_Gradient中，问题是什么意义
			mod = *(pt + xs * xo + ys * yo);		//显然这是梯度，模长越长，落差越大
			ang = *(pt + xs * xo + ys * yo + 1);	//这个是角度，该点的梯度方向
			fbin = (float)(nbins * ang / (2 * PI));			//一个圆分36分，fbin即ang在哪一份

#if defined(VL_SIFT_BILINEAR_ORIENTATIONS)
			{
				int bin = (int)floor(fbin - 0.5f);
				float rbin = fbin - bin - 0.5f;
				//一个点的角度一般跨了两份角度，按照位置决定哪份多点哪份少点，然后
				//梯度的贡献也分开两部分，分配到相邻的两个角度上。这个计算就是个小学算数
				//画个图，如果落在两份之间，当然是各占0.5。如果向右多点，则左邻居少点。
				//按比例一分完事

				// 此处傻逼，摆了个乌龙，存在负位置
				//hist[(bin + nbins) % nbins] += (1 - rbin) * mod * wgt;
				//hist[(bin + 1) % nbins] += (rbin)*mod * wgt;

				hist[(bin + nbins) % nbins] += (1 - rbin) * mod * wgt;
				hist[(bin + 1 + nbins) % nbins] += (rbin)*mod * wgt;

			}
#else
			{
				int    bin = vl_floor_d(fbin);
				bin = vl_floor_d(nbins * ang / (2 * PI));
				hist[(bin) % nbins] += mod * wgt;
			}
#endif

		} /* for xs */
	} /* for ys */

	/* smooth histogram */
	for (iter = 0; iter < 6; iter++)
	{
		float prev = hist[nbins - 1];
		float first = hist[0];
		int i;
		//printf("iter:%d\n", iter);
		for (i = 0; i < nbins - 1; i++)
		{
			//以自己为中心，加上左右邻居，再来个平均作为新值
			//然后来回折腾几次已达到平滑效果
			float newh = (prev + hist[i] + hist[(i + 1) % nbins]) / 3.0f;
			prev = hist[i];
			hist[i] = newh;
		}
		hist[i] = (prev + hist[i] + first) / 3.0f;
	}

	/* find the histogram maximum */
	maxh = 0;
	for (i = 0; i < nbins; ++i)
		maxh = max(maxh, hist[i]);

	/* find peaks within 80% from max */
	nangles = 0;
	for (i = 0; i < nbins; ++i)
	{
		float h0 = hist[i];
		float hm = hist[(i - 1 + nbins) % nbins];	//h0角度旁边的角度
		float hp = hist[(i + 1 + nbins) % nbins];

		/* is this a peak? */
		if (h0 > 0.8 * maxh && h0 > hm && h0 > hp)
		{//由于W足够大，表示关于本关键点的搜索范围足够大，比如半径达到8像素，这个时候完全可以右多个落差很大的路径
			//想都想到，下山的路径不止一条。那么这个时候接着看梯度是否足够大，如果大到一定程度，这里是最大梯度的0.8倍，
			//就能形成另一各特征,这个特征也是用方向来表达
			/* quadratic interpolation */
			float di = -0.5f * (hp - hm) / (hp + hm - 2.f * h0);
			float th = (float)(2.f * PI * (i + di + 0.5f) / nbins);
			Angle[nangles++] = th;
			if (nangles == 4)
				goto enough_angles;
		}
	}
enough_angles:
	return nangles;
}

inline float Normalize_Histogram(float* begin, float* end)
{//基本上就是向量的单位化
#define eps 1.19209290E-07F
	float* iter;
	float  norm = 0.0;

	union {
		__m512 Total_16;
		__m256 Total_8[2];
		__m128 Total_4[2];
	};
	__m512 Norm_Recip;
	Total_16 = _mm512_mul_ps(*(__m512*)begin, *(__m512*)begin);
	for (iter = begin + 16; iter < end; iter += 16)
		Total_16 = _mm512_add_ps(Total_16, _mm512_mul_ps(*(__m512*)iter, *(__m512*)iter));
	Total_8[0] = _mm256_add_ps(Total_8[0], Total_8[1]);
	Total_4[0] = _mm_add_ps(Total_4[0], Total_4[1]);
	//norm = Total_4[0].m128_f32[0] + Total_4[0].m128_f32[1] + Total_4[0].m128_f32[2] + Total_4[0].m128_f32[3];
	norm = ((float*)&Total_4[0])[0] + ((float*)&Total_4[0])[1] + ((float*)&Total_4[0])[2] + ((float*)&Total_4[0])[3];
	norm = (float)sqrt(norm);
	Norm_Recip = _mm512_set1_ps(1.f / norm);

	for (iter = begin; iter != end; iter += 16)
		*(__m512*)iter = _mm512_mul_ps(*(__m512*)iter, Norm_Recip);
	return norm;
#undef eps
}


void Calc_Keypoint_Descriptor_AVX512(Sift* poSift, float* descr, Sift_Keypoint const* k, float angle0)
{//简单加速
#define eps 2.220446049250313e-16
#define NBO 8
#define NBP 4
	float const magnif = poSift->magnif;
	float       xper = (float)pow(2.0f, poSift->o_cur);
	int          w = poSift->octave_width;
	int          h = poSift->octave_height;
	int const    xo = 2;         /* x-stride */
	int const    yo = 2 * w;     /* y-stride */
	//int const    so = 2 * w * h; /* s-stride */
	float       x = k->x / xper;
	float       y = k->y / xper;
	float       sigma = k->sigma / xper;

	int          xi = (int)(x + 0.5);	//此处xi用的是子像素值，但是又四舍五入为整数，似乎没什么用
	int          yi = (int)(y + 0.5);
	int          si = k->is;

	float const st0 = (float)sin(angle0);
	float const ct0 = (float)cos(angle0);
	float const SBP = (float)(magnif * sigma + eps);
	//以下这个W就是r
	int    const W = (int)(floor(sqrt(2.0) * SBP * (NBP + 1) / 2.0 + 0.5));

	int const binto = 1;          /* bin theta-stride */
	int const binyo = NBO * NBP;  /* bin y-stride */
	int const binxo = NBO;        /* bin x-stride */

	int bin, dxi, dyi;
	float* pt;
	float* dpt;

	//越界检查
	if (k->o != poSift->o_cur || xi < 0 || xi >= w ||
		yi < 0 || yi >= h - 1 || si    <  poSift->s_min + 1 ||
		si    >  poSift->s_max - 2)
		return;

	//8x4x4块，128维Descriptor，要进一步分析其划分
	memset(descr, 0, sizeof(float) * NBO * NBP * NBP);

	//x0：每点耗费多少个样本， y0：每行多少个样本。 因为每点的梯度由两个样本构成
	//1，落差；2，方向  pt的含义是：point
	pt = &poSift->m_pGrad[xi * xo + yi * yo];

	//从以下这段猜测，一个Descriptor和 Surf一样，右4x4个sub block构成，每个
	//sub block的的大小为2x2=4, 所以一共有4x4lx4组信息。每组信息可能由两个字段
	//构成： 1，模长（落差）；2，方向。 所以一共128维
	//dpt的含义是 descriptor point，从数值上，dpt是位于descr[80], 为什么是80？
	//这是因为这个位置是在4x4的中间块， [2][2]上。 2*4*2*2 + 2*2*2=80，齐活
	dpt = descr + (NBP / 2) * binyo + (NBP / 2) * binxo;

#define atd(dbinx,dbiny,dbint) *(dpt + (dbint)*binto + (dbiny)*binyo + (dbinx)*binxo)

	//由于是简单优化，前面照抄，后面并行操作再优化
	int y_Start = max(-W, 1 - yi),
		y_End = min(+W, h - yi - 2),
		x_Start = max(-W, 1 - xi),
		x_End = min(+W, w - xi - 2),//iWin_Width = x_End - x_Start + 1,
		y1;
	float* pt_Cur, dy,
		SBP_recip = 1.f / SBP,
		st0_x_dy,	//st0*dy
		ct0_x_dy,	//ct0*dy
		wsigma_recip = 1.f / (2.0f * poSift->windowSize * poSift->windowSize);

	int iResult, dbinx, dbiny, dbint, iRemain;
	int iResult_x, iResult_y;
	__m512 mod_16, dx_16, nx_16, ny_16, win_16,
		rbinx_16, rbiny_16, rbint_16,
		weight_16;
	union { __m512 angle_16, theta_16, nt_16; };
	__m512i Index_16, dxi_16, binx_dbinx, biny_dbiny, bint_dbint,
		binx_16, biny_16, bint_16,
		seq_16 = _mm512_set_epi32(15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0);

	//整个方块窗口尽可能大
	for (dyi = y_Start; dyi <= y_End; dyi++)
	{
		pt_Cur = pt + dyi * yo;
		y1 = dyi * yo;
		//dy_16 = _mm512_set1_ps(yi + dyi - y);
		dy = yi + dyi - y;
		st0_x_dy = st0 * dy;
		ct0_x_dy = ct0 * dy;
		for (dxi = x_Start; dxi <= x_End; dxi += 16, pt_Cur += 16)
		{
			iRemain = (x_End - dxi + 1);
			if (iRemain > 16)
				iRemain = 16;

			//float mod = *(pt + dxi * xo + dyi * yo + 0);
			//Index_16 = _mm512_add_epi32(dxi_512, _mm512_set1_epi32(dxi));
			//Index_16 = _mm512_mullo_epi32(Index_16, _mm512_set1_epi32(xo));
			//Index_16 = _mm512_add_epi32(_mm512_set1_epi32(y1), Index_16);
			dxi_16 = _mm512_add_epi32(seq_16, _mm512_set1_epi32(dxi));
			Index_16 = _mm512_add_epi32(_mm512_set1_epi32(y1), _mm512_mullo_epi32(dxi_16, _mm512_set1_epi32(xo)));
			mod_16 = _mm512_i32gather_ps(Index_16, pt, 4);

			//Index_16=_mm512_add_epi32(Index_16, _mm512_set1_epi32(1));
			angle_16 = _mm512_i32gather_ps(_mm512_add_epi32(Index_16, _mm512_set1_epi32(1)), pt, 4);

			//float theta = vl_mod_2pi_f(angle - angle0);
			angle_16 = _mm512_sub_ps(angle_16, _mm512_set1_ps(angle0));
			iResult = _mm512_cmp_ps_mask(angle_16, _mm512_set1_ps(0), _CMP_LT_OS);
			theta_16 = _mm512_mask_add_ps(angle_16, iResult, angle_16, _mm512_set1_ps((float)(PI * 2)));

			//float dx = (float)(xi + dxi - x);
			dx_16 = _mm512_sub_ps(_mm512_cvtepi32_ps(_mm512_add_epi32(dxi_16, _mm512_set1_epi32(xi))), _mm512_set1_ps(x));

			//float nx = (ct0 * dx + st0 * dy) / SBP;
			//nx_16 = _mm512_mul_ps(dx_16, _mm512_set1_ps(ct0));
			//nx_16 = _mm512_add_ps(nx_16, st0_x_dy);
			//nx_16 = _mm512_mul_ps(nx_16, _mm512_set1_ps(SBP_recip));
			nx_16 = _mm512_mul_ps(_mm512_add_ps(_mm512_mul_ps(dx_16, _mm512_set1_ps(ct0)), _mm512_set1_ps(st0_x_dy)), _mm512_set1_ps(SBP_recip));

			//float ny = (-st0 * dx + ct0 * dy) / SBP;
			//ny_16 = _mm512_mul_ps(dx_16, _mm512_set1_ps(-st0));
			//ny_16 = _mm512_add_ps(ny_16, _mm512_set1_ps(ct0_x_dy));
			//ny_16 = _mm512_mul_ps(ny_16, _mm512_set1_ps(SBP_recip));
			ny_16 = _mm512_mul_ps(_mm512_add_ps(_mm512_mul_ps(dx_16, _mm512_set1_ps(-st0)), _mm512_set1_ps(ct0_x_dy)), _mm512_set1_ps(SBP_recip));

			//float nt = NBO * theta / (2 * PI);
			//nt_16 = _mm512_mul_ps(_mm512_set1_ps(NBO), theta_16);
			//nt_16 = _mm512_mul_ps(nt_16, _mm512_set1_ps(1.f / (2.f * PI)));
			nt_16 = _mm512_mul_ps(_mm512_mul_ps(_mm512_set1_ps(NBO), theta_16), _mm512_set1_ps((float)(1.f / (2.f * PI))));

			//float win = fast_expn((nx * nx + ny * ny) / (2.0f * wsigma * wsigma));
			win_16 = _mm512_mul_ps(_mm512_add_ps(_mm512_mul_ps(nx_16, nx_16), _mm512_mul_ps(ny_16, ny_16)), _mm512_set1_ps(-wsigma_recip));
			win_16 = _mm512_exp_ps(win_16);

			//int binx = (int)vl_floor_f(nx - 0.5f);
			//int         biny = (int)vl_floor_f(ny - 0.5f);
			//int         bint = (int)vl_floor_f(nt);
			binx_16 = _mm512_cvtps_epi32(_mm512_floor_ps(_mm512_sub_ps(nx_16, _mm512_set1_ps(0.5f))));
			biny_16 = _mm512_cvtps_epi32(_mm512_floor_ps(_mm512_sub_ps(ny_16, _mm512_set1_ps(0.5f))));
			bint_16 = _mm512_cvtps_epi32(_mm512_floor_ps(nt_16));

			//float rbinx = nx - (binx + 0.5f);
			//float rbiny = ny - (biny + 0.5f);
			//float rbint = nt - bint;
			rbinx_16 = _mm512_sub_ps(nx_16, _mm512_add_ps(_mm512_cvtepi32_ps(binx_16), _mm512_set1_ps(0.5f)));
			rbiny_16 = _mm512_sub_ps(ny_16, _mm512_add_ps(_mm512_cvtepi32_ps(biny_16), _mm512_set1_ps(0.5f)));
			rbint_16 = _mm512_sub_ps(nt_16, _mm512_cvtepi32_ps(bint_16));
			for (dbinx = 0; dbinx < 2; ++dbinx)
			{
				//if (!(binx + dbinx >= -(NBP / 2) && binx + dbinx < (NBP / 2)))
					//continue;
				binx_dbinx = _mm512_add_epi32(binx_16, _mm512_set1_epi32(dbinx));
				iResult_x = _mm512_cmp_epi32_mask(binx_dbinx, _mm512_set1_epi32(-(NBP / 2)), _MM_CMPINT_GE);
				iResult_x &= _mm512_cmp_epi32_mask(binx_dbinx, _mm512_set1_epi32(NBP / 2), _MM_CMPINT_LT);
				if (!iResult_x)
					continue;

				for (dbiny = 0; dbiny < 2; ++dbiny)
				{
					//if (!(biny + dbiny >= -(NBP / 2) && biny + dbiny < (NBP / 2)))
						//continue;
					biny_dbiny = _mm512_add_epi32(biny_16, _mm512_set1_epi32(dbiny));
					iResult_y = iResult_x & _mm512_cmp_epi32_mask(biny_dbiny, _mm512_set1_epi32(-(NBP / 2)), _MM_CMPINT_GE);
					iResult_y &= _mm512_cmp_epi32_mask(biny_dbiny, _mm512_set1_epi32(NBP / 2), _MM_CMPINT_LT);
					if (!iResult_y)
						continue;
					for (dbint = 0; dbint < 2; dbint++)
					{
						iResult = (iResult_x & iResult_y) & ((1 << iRemain) - 1);
						//一条指令顶半边天
						//bint_dbint = _mm512_rem_epi32(_mm512_add_epi32(bint_16, _mm512_set1_epi32(dbint)), _mm512_set1_epi32(NBO));
						bint_dbint = _mm512_and_epi32(_mm512_add_epi32(bint_16, _mm512_set1_epi32(dbint)), _mm512_set1_epi32(0x7));
						weight_16 = _mm512_abs_ps(_mm512_sub_ps(_mm512_set1_ps(1.f - dbinx), rbinx_16));
						weight_16 = _mm512_mul_ps(weight_16, _mm512_abs_ps(_mm512_sub_ps(_mm512_set1_ps(1.f - dbiny), rbiny_16)));
						weight_16 = _mm512_mul_ps(weight_16, _mm512_abs_ps(_mm512_sub_ps(_mm512_set1_ps(1.f - dbint), rbint_16)));
						weight_16 = _mm512_mul_ps(weight_16, _mm512_mul_ps(win_16, mod_16));

						Index_16 = _mm512_mullo_epi32(binx_dbinx, _mm512_set1_epi32(binxo));
						Index_16 = _mm512_add_epi32(Index_16, _mm512_mullo_epi32(biny_dbiny, _mm512_set1_epi32(binyo)));
						Index_16 = _mm512_add_epi32(Index_16, _mm512_mullo_epi32(bint_dbint, _mm512_set1_epi32(binto)));

						//iResult &= (1 << iRemain) - 1;						
						for (int i = 0; iResult;)
						{
							if (iResult & 1)
								*(dpt + ((int*)&Index_16)[i]) += ((float*)&weight_16)[i];
							//*(dpt + Index_16.m512i_i32[i]) += weight_16.m512_f32[i];

							iResult >>= 1;
							i++;
						}
					}
				}
			}
		}
	}

	//Disp(descr, 128, 1);
	float norm = Normalize_Histogram(descr, descr + NBO * NBP * NBP);
	/* Set the descriptor to zero if it is lower than our norm_threshold */
	if (poSift->norm_thresh && norm < poSift->norm_thresh)
	{
		for (bin = 0; bin < NBO * NBP * NBP; ++bin)
			descr[bin] = 0;
	}
	else
	{
		for (bin = 0; bin < NBO * NBP * NBP; bin += 16)
		{
			iResult = _mm512_cmp_ps_mask(*(__m512*)&descr[bin], _mm512_set1_ps(0.2f), _MM_CMPINT_GE);
			if (iResult)
				_mm512_mask_i32scatter_ps(&descr[bin], iResult, seq_16, _mm512_set1_ps(0.2f), 4);
		}
		Normalize_Histogram(descr, descr + NBO * NBP * NBP);
	}
	return;
#undef eps
#undef NBO
#undef NBP
}

void L1RootNormalizeFeatureDescriptors(float Desc[])
{//对128维Descriptor进行另一种规格化。 模为一维范数
	float fFactor;
	int i;
	for (i = 0, fFactor = 0; i < 128; i++)
		fFactor += Desc[i];
	fFactor = 1.f / fFactor;
	for (i = 0; i < 128; i++)
		Desc[i] = (float)sqrt(Desc[i] * fFactor);
}

void FeatureDescriptorsToUnsignedByte(float Desc[128], unsigned char Desc_i[128])
{//只求数字一致，不管别的	0, 7, 6, 5, 4, 3, 2, 1, 8, 15, 14, 13, 12, 11, 10, 9
	__m512i Re_Order = _mm512_set_epi32(9, 10, 11, 12, 13, 14, 15, 8, 1, 2, 3, 4, 5, 6, 7, 0);
	__m512 Value_16 = _mm512_set1_ps(512.f);

	for (int i = 0; i < 128; i += 16)
		*(__m128i*)& Desc_i[i] = _mm512_cvtepi32_epi8(_mm512_cvt_roundps_epi32(_mm512_mul_ps(_mm512_i32gather_ps(Re_Order, &Desc[i], 4), Value_16), _MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC));
}

static int Get_Feature(Sift* poSift, Sift_Feature_Level* poFeature_Level, int iDoG_Top, int iOctave_Top)
{//继续向上走
	union {
		struct {
			float* pt, * pOctave_0, * pOctave_1;
		};
		__m256i Addr_3;
	};
	Sift oSift = *poSift;
	Sift_Feature_Level oFeature_Level = *poFeature_Level;
	Sift_Keypoint oKey_Point;
	Sift_Feature* poFeature_Byte;
	//Sift_Feature oFeature;
	float Desc[128];
	float* pEnd, * pFilter, angles[4];;
	float sd;
	int i, r, num_orientations, num_used_orientations, iRemain,
		iSize = oSift.octave_width * oSift.octave_height;
	int iTemp, iOctave_Top_Index = iOctave_Top % 3;	//根据当前处理的最高层来推出在内存数组中的索引
	int bRet = 0;

	//做一层高斯滤波
	pOctave_0 = oSift.m_pGauss[iOctave_Top_Index];
	iTemp = (1 + iOctave_Top_Index) % 3;
	pOctave_1 = oSift.m_pGauss[iTemp];
	oSift.m_pTemp = oSift.m_pDoG[iOctave_Top_Index];	//共用本层DoG
	sd = oSift.dsigma0 * (float)pow(oSift.sigmak, iOctave_Top);
	r = (int)max((float)ceil(4.0 * sd), 1.f);
	Gen_Gauss_Filter(r, sd, &pFilter);
	if (!pFilter)
		goto END;

	//从下面一层生成上面一层
	Gauss_Filter_AVX512(pOctave_0, oSift.octave_width, oSift.octave_height, pOctave_1, oSift.m_pTemp, r, pFilter);
	Free(pFilter);

	//做一层DoG
	pt = oSift.m_pDoG[(iDoG_Top + 1) % 3];
	pEnd = pt + iSize;
	for (; pt < pEnd; Addr_3 = _mm256_add_epi64(Addr_3, _mm256_set1_epi64x(16 * 4))   /*pt += 16, pOctave_0 += 16, pOctave_1 += 16*/)
		*(__m512*)pt = _mm512_sub_ps(*(__m512*)pOctave_1, *(__m512*)pOctave_0);

	//对已经做好的DoG做一次Keypoint检测
	Get_Keypoint(&oSift, iDoG_Top, iOctave_Top - 2, &oFeature_Level);
	if (iOctave_Top == oSift.s_max)
	{//若当前的高斯层已经到了最高层，则将当前层做一个下采样，存在下一层
		Resize_Image_Half(oSift.m_pGauss[(iOctave_Top - 1) % 3], oSift.octave_width, oSift.octave_height, oSift.m_pGauss[iOctave_Top % 3]);
		//bSave_Image("c:\\tmp\\1.bmp", oSift.m_pGauss[iOctave_Top % 3], oSift.octave_width>>1, oSift.octave_height>>1);
	}
	//更新梯度, 其实应该用Octave[iOctave_Top_Index-1]这帧，由于减法-1 再求模会有负数问题，
	//所以将 -1求模改成 +2 %3，对等
	pOctave_0 = oSift.m_pGauss[iDoG_Top % 3];
	oSift.m_pGrad = oSift.m_pDoG[(iDoG_Top - 1) % 3];
	_Update_Gradiant_AVX512(pOctave_0, oSift.m_pGrad, oSift.octave_width, oSift.octave_height);

	static int iCounter = 0;
	iRemain = oSift.m_iMax_Feature_Count - (int)(oFeature_Level.m_pFeature - oSift.m_pFeature);	//还剩下读诵好个
	for (i = 0; i < oFeature_Level.m_iKeypoint_Count; i++)
	{
		oKey_Point = oFeature_Level.m_pKeypoint[i];
		//if (iCounter >= 235)
			//printf("here\n");
		num_orientations = iCal_Key_Point_Oriatation(&oSift, angles, oKey_Point);
		num_used_orientations = min(num_orientations, 2);
		for (int o = 0; o < num_used_orientations; o++)
		{
			if (iRemain <= 0)
			{//需要扩大Size
				if (!(Resize_Keypoint_Feature_Size(&oSift, &oFeature_Level, 0)))
					return 0;
				*poSift = oSift;
				iRemain += oSift.m_iMax_Keypoint_Count >> 1;
			}
			iRemain--;

			poFeature_Byte = &oFeature_Level.m_pFeature[oFeature_Level.m_iFeature_Count];
			poFeature_Byte->x = oKey_Point.x + 0.5f;
			poFeature_Byte->y = oKey_Point.y + 0.5f;
			poFeature_Byte->x_Keypoint = oKey_Point.ix;
			poFeature_Byte->y_Keypoint = oKey_Point.iy;
			poFeature_Byte->angle = angles[o];
			poFeature_Byte->o = oSift.o_cur;;
			poFeature_Byte->s = oKey_Point.is;
			Calc_Keypoint_Descriptor_AVX512(&oSift, Desc, &oKey_Point, angles[o]);
			L1RootNormalizeFeatureDescriptors(Desc);
			FeatureDescriptorsToUnsignedByte(Desc, poFeature_Byte->m_Desc_i);
			//Disp_Feature(*poFeature_Byte);
			
			//if (abs(poFeature_Byte->x - 4.439735e+02f) < 0.0001f)
				//Disp(poFeature_Byte->m_Desc_i, 1, 128, "Feature");

				//printf("here");
			//printf("%e\n", poFeature_Byte->x);
			iCounter++;
			oFeature_Level.m_iFeature_Count++;
		}
	}
	*poFeature_Level = oFeature_Level;
	*poSift = oSift;
	bRet = 1;
END:
	return bRet;
}

static void Update_Sift(Sift* poSift, int iOctave_Reserve_Index)
{//此时要将Octave[iOctave_Reserve_Index]放在Octave[0]处，并缩小一半
	Sift oSift = *poSift;
	//Resize_Image_Half(oSift.m_pGauss[iOctave_Reserve_Index], oSift.octave_width, oSift.octave_height, oSift.m_pGauss[0]);
	oSift.o_cur += 1;
	oSift.octave_width >>= 1, oSift.octave_height >>= 1;
	Arrange_Mem(&oSift, oSift.octave_width, oSift.octave_height, oSift.m_pBuffer_Start);
	*poSift = oSift;
	return;
}

int Get_Sift_Feature(Image oImage, 	float (**ppKey_Point)[2], unsigned char (**ppFeature)[128], int* piFeature_Count, int o_min, float fPeak_Thresh)
{//Sift计算接口，尽可能清爽。 返回特征集及数量
	Sift oSift;
	Sift_Feature_Level oFeature_Level = {}, Feature_Level[4][3];
	int i, iFeature_Count = 0, iKeypoint_Count = 0,bRet=1;

	Init_Sift(&oSift, 4, 3, o_min, oImage.m_iWidth, oImage.m_iHeight, fPeak_Thresh);
	if (!oSift.m_pBuffer_Start)
		goto END;

	Get_w_h(oSift.m_iWidth, oSift.m_iHeight, oSift.octave_width, oSift.octave_height, oSift.o_cur);
	Arrange_Mem(&oSift, oSift.octave_width, oSift.octave_height, oSift.m_pBuffer_Start);
	oSift.m_pKey_Point = oFeature_Level.m_pKeypoint = (Sift_Keypoint*)pMalloc(oSift.m_iMax_Feature_Count * sizeof(Sift_Keypoint));
	oSift.m_pFeature = oFeature_Level.m_pFeature = (Sift_Feature*)pMalloc(oSift.m_iMax_Feature_Count * sizeof(Sift_Feature));
	if (!oFeature_Level.m_pKeypoint || !oFeature_Level.m_pFeature)
		goto END;

	float* pImage;
	if (oSift.o_min == 0)
	{
		pImage = oSift.m_pGauss[0];
		RGB_2_Gray(oImage, pImage);
	}else
	{
		pImage = (float*)pMalloc(oImage.m_iWidth * oImage.m_iHeight * sizeof(float));
		RGB_2_Gray(oImage, pImage);
		if (oSift.o_min == -1)
			Resize_Image_x2(pImage, oSift.m_iWidth, oSift.m_iHeight, oSift.m_pGauss[0]);
		else if (oSift.o_min == 1)
			Resize_Image_Half(pImage, oSift.m_iWidth, oSift.m_iHeight, oSift.m_pGauss[0]);
		Free(pImage);
	}

	while (1)
	{//逐个分辨率做Sift分析，w*2*h*2->w*h->w/2 * h/2 -> w/4 * h/4
		//此处要一改前人方法，省点内存
		Sift_Dectact(&oSift, 0);
		for (i = 0; i < 3; i++)
		{
			Get_Feature(&oSift, &oFeature_Level, i + 1, i + 2);

			Feature_Level[oSift.o_cur + 1][i] = oFeature_Level;
			oFeature_Level.m_pFeature += oFeature_Level.m_iFeature_Count;
			oFeature_Level.m_pKeypoint += oFeature_Level.m_iKeypoint_Count;
			//printf("o:%d s:%d Keypoint:%d Feature:%d\n", oSift.o_cur, i, oFeature_Level.m_iKeypoint_Count,oFeature_Level.m_iFeature_Count);
			iFeature_Count += oFeature_Level.m_iFeature_Count;
			iKeypoint_Count += oFeature_Level.m_iKeypoint_Count;
			oFeature_Level.m_iFeature_Count = oFeature_Level.m_iKeypoint_Count = 0;
		}
		if (oSift.o_cur == 2)	//目前只看到最多搞到第2层octave，还要进一步验证
			break;
		//更新一下oSift
		Update_Sift(&oSift, 0);
	}
	if (piFeature_Count)
		*piFeature_Count = iFeature_Count;

	if (ppFeature)
	{
		unsigned char (*pFeature)[128] = (unsigned char(*)[128])pMalloc(iFeature_Count * 128);
		for (int i = 0; i < iFeature_Count; i++)
			memcpy(pFeature[i], oSift.m_pFeature[i].m_Desc_i, 128);
		*ppFeature = pFeature;
	}

	if (ppKey_Point)
	{
		float (*pKey_Point)[2] = (float(*)[2])pMalloc(iFeature_Count * 2 * sizeof(float));
		for (int i = 0; i < iFeature_Count; i++)
		{
			pKey_Point[i][0] = oSift.m_pFeature[i].x;	// oSift.m_pKey_Point[i].x;
			pKey_Point[i][1] = oSift.m_pFeature[i].y;	// oSift.m_pKey_Point[i].y;
			if (abs(pKey_Point[i][0] - 4.439735e+02f) < 0.0001f)
				printf("%f\n", pKey_Point[i][0]);
		}
		*ppKey_Point = pKey_Point;
	}
	bRet = 1;
END:
	Free(oSift.m_pBuffer_Start);
	Free(oSift.m_pKey_Point);
	Free(oSift.m_pFeature);
	if (bRet)
	{//将数据提前，尽可能少点碎片
		if (ppFeature)
			*ppFeature = (unsigned char (*)[128])pMove_2_Front(*ppFeature);
		if (ppKey_Point)
			*ppKey_Point = (float(*)[2]) pMove_2_Front(*ppKey_Point);
	}
	return 1;
}

int Get_Sift_Feature(Image oImage,  Point_2D<float>** ppKey_Point,
	unsigned char (**ppFeature)[128], int* piFeature_Count, int o_min, float fPeak_Thresh)
{
	int iResult = 0, iCount;	
	if(ppKey_Point)
	{
		float (*pKey_Point)[2];
		iResult = Get_Sift_Feature(oImage, &pKey_Point, ppFeature, &iCount, o_min, fPeak_Thresh);
		Point_2D<float>* pKey_Point_1 = (Point_2D<float>*)pMalloc(iCount * sizeof(Point_2D<float>));
		for (int i = 0; i < iCount; i++)
		{
			pKey_Point_1[i].m_iCamera_Index = -1;
			pKey_Point_1[i].m_iPoint_Index = i;
			pKey_Point_1[i].m_Pos[0] = pKey_Point[i][0];
			pKey_Point_1[i].m_Pos[1] = pKey_Point[i][1];
		}
		Free(pKey_Point);
		*ppKey_Point = (Point_2D<float>*)pMove_2_Front(pKey_Point_1);;
	}else
		iResult = Get_Sift_Feature(oImage, (float(**)[2])NULL, ppFeature, &iCount, o_min);
	*piFeature_Count = iCount;
	return iResult;
}

int Get_Sift_Feature_1(const char* pcFile, float(**ppKey_Point)[2], unsigned char (**ppFeature)[128], int* piCount, int o_min, float fPeak_Thresh)
{//Sift 获得
	Image oImage;	
	if (!bLoad_Image(pcFile, &oImage))
		return 0;
	int iResult = Get_Sift_Feature(oImage, ppKey_Point,ppFeature, piCount, o_min, fPeak_Thresh);
	Free_Image(&oImage);
	//Disp_Mem();
	return iResult;
}

int Get_Sift_Feature(const char* pcFile, Point_2D<float> **ppKey_Point, unsigned char (**ppFeature)[128], int* piCount, int o_min, float fPeak_Thresh)
{//Sift 获得
	Image oImage;
	if (!bLoad_Image(pcFile, &oImage))
		return 0;
	unsigned long long tStart = iGet_Tick_Count();
	int iResult = Get_Sift_Feature(oImage,ppKey_Point, ppFeature, piCount, o_min, fPeak_Thresh);
	Free_Image(&oImage);
	printf("%lld\n", iGet_Tick_Count() - tStart);
	return iResult;
}
template<typename _T>Match_Pose* Pose_Base_Get_Match_Poses(Pose_Base<_T>* poBase, int iPose_ID,int *piCount )
{//给定一个Pose_ID，返回一组与此Pose匹配的Pose
	Pose_Bsse_Pose<_T>* poPose = &poBase->m_pPose[iPose_ID];
	*piCount = poPose->m_oMatch_Pose_Start_Count.m_iCount;
	if (!poPose->m_oMatch_Pose_Start_Count.m_iCount)
		return NULL;
	return &poBase->m_pMatch_Pose[poPose->m_oMatch_Pose_Start_Count.m_iStart];
}

template<typename _T>Match_Pose* Pose_Base_Get_Match_Pose(Pose_Base<_T>*poBase, int iPose_A, int iPose_B)
{
	int iCount;
	Match_Pose* pPose = Pose_Base_Get_Match_Poses(poBase, iPose_A, &iCount);
	if (!pPose)
		return NULL;
	for (int i = 0; i < iCount; i++)
	{
		if (pPose[i].m_iPose_B == iPose_B)
			return &pPose[i];
	}
	return NULL;
}
template int Pose_Base_Add_Match(Pose_Base<double>* poBase, int iPose_A, int iPose_B, unsigned short Match_Index[][2], unsigned short int iMatch_Index_Count);
template<typename _T>int Pose_Base_Add_Match(Pose_Base<_T>* poBase, int iPose_A,int iPose_B, unsigned short Match_Index[][2], unsigned short int iMatch_Index_Count)
{
	const int iFill_Count = 2;
	int bSwap_Index = 0;
	if (iPose_A > iPose_B)
	{
		SWAP(int, iPose_A, iPose_B);
		bSwap_Index = 1;
	}

	Start_Count *poA_Match_Pose_Start_Count = &poBase->m_pPose[iPose_A].m_oMatch_Pose_Start_Count;
	Start_Count *poB_Match_Pose_Start_Count = &poBase->m_pPose[iPose_B].m_oMatch_Pose_Start_Count;
	//保守一点，先确顶是否有足够内存，精确计算太麻烦，大体是全部新加入是否足够
	{
		if (poBase->m_iMatch_Pose_Cur + poA_Match_Pose_Start_Count->m_iCount + iFill_Count +
			poB_Match_Pose_Start_Count->m_iCount + iFill_Count >= poBase->m_iMax_Match_Pose_Count)
		{//未雨绸缪，直接加内存
			int iMax_Match_Pose_Count = max(poBase->m_iMatch_Pose_Cur + poA_Match_Pose_Start_Count->m_iCount + 1 +
				poB_Match_Pose_Start_Count->m_iCount + 1, poBase->m_iMax_Match_Pose_Count * 2);
			if (!bExpand_Ex((void**)&poBase->m_pMatch_Pose, iMax_Match_Pose_Count * sizeof(Match_Pose)))
				return 0;
			poBase->m_iMax_Match_Pose_Count = iMax_Match_Pose_Count;
		}
		if (poBase->m_iMatch_Index_Cur + iMatch_Index_Count >= poBase->m_iMax_Match_Index_Count)
		{
			int iMax_Match_Index_Count = max(poBase->m_iMatch_Index_Cur + iMatch_Index_Count, poBase->m_iMax_Match_Index_Count * 2);
			if (!bExpand_Ex((void**)&poBase->m_pMatch_Index, iMax_Match_Index_Count * 2 * sizeof(unsigned short)))
				return 0;
			poBase->m_iMax_Match_Index_Count = iMax_Match_Index_Count;
		}
	}

	//先确定是否重合
	int iCount;
	Match_Pose* poCur, * poA_Match_Poses, * poB_Match_Poses;
	poCur = poA_Match_Poses = Pose_Base_Get_Match_Poses(poBase, iPose_A, &iCount);
	for (int i = 0; i < iCount; i++, poCur++)
	{
		if (poCur->m_iPose_B == iPose_B)
		{
			printf("Match Pose %d and %d existed\n", iPose_A, iPose_B);
			return 0;
		}
	}
	
	//将match_Index 加到最后
	unsigned short (*pCur)[2] = &poBase->m_pMatch_Index[poBase->m_iMatch_Index_Cur];
	Start_Count oMatch_Index_Start_Count = {(int)poBase->m_iMatch_Index_Cur ,iMatch_Index_Count };
	if (bSwap_Index)
		for (int i = 0; i < iMatch_Index_Count; i++)
			pCur[i][0] = Match_Index[i][1],
			pCur[i][1] = Match_Index[i][0];
	else
		memcpy(pCur, Match_Index, iMatch_Index_Count * 2 * sizeof(unsigned short));
	poBase->m_iMatch_Index_Cur += iMatch_Index_Count;


	//看看能否塞入桶中
	Match_Pose* pNew;
	if (!poA_Match_Poses)
	{
		*poA_Match_Pose_Start_Count = { (int)poBase->m_iMatch_Pose_Cur,1 };
		pNew = &poBase->m_pMatch_Pose[poBase->m_iMatch_Pose_Cur];
		poBase->m_iMatch_Pose_Cur += iFill_Count;
	}
	else if (poA_Match_Poses[iCount].m_oMatch_Index_Start_Count.m_iCount)
	{//已经占，只能异地加入
		*poA_Match_Pose_Start_Count = { (int)poBase->m_iMatch_Pose_Cur,poA_Match_Pose_Start_Count->m_iCount + 1 };
		pNew = &poBase->m_pMatch_Pose[poBase->m_iMatch_Pose_Cur];
		memcpy(pNew, poA_Match_Poses, iCount * sizeof(Match_Pose));
		poBase->m_iMatch_Pose_Cur += iCount + iFill_Count;
		pNew += iCount;
	}
	else		//可以直接加入
	{
		poA_Match_Pose_Start_Count->m_iCount++;
		pNew = &poBase->m_pMatch_Pose[poA_Match_Pose_Start_Count->m_iStart + poA_Match_Pose_Start_Count->m_iCount];
	}

	pNew->m_iPose_B = iPose_B;
	pNew->m_oMatch_Index_Start_Count = oMatch_Index_Start_Count;

	//再轮到Pose_B
	poCur = poB_Match_Poses = Pose_Base_Get_Match_Poses(poBase, iPose_A, &iCount);
	if (!poB_Match_Poses)
	{
		*poB_Match_Pose_Start_Count = { (int)poBase->m_iMatch_Pose_Cur,1 };
		pNew = &poBase->m_pMatch_Pose[poBase->m_iMatch_Pose_Cur];
		poBase->m_iMatch_Pose_Cur += iFill_Count;
	}else if (poB_Match_Poses[iCount].m_oMatch_Index_Start_Count.m_iCount)
	{//已经占，只能异地加入
		*poB_Match_Pose_Start_Count = { (int)poBase->m_iMatch_Pose_Cur,poB_Match_Pose_Start_Count->m_iCount + 1 };
		pNew = &poBase->m_pMatch_Pose[poBase->m_iMatch_Pose_Cur];
		memcpy(pNew, poB_Match_Poses, iCount * sizeof(Match_Pose));
		poBase->m_iMatch_Pose_Cur += iCount + iFill_Count;
		pNew += iCount;
	}else		//可以直接加入
	{
		poB_Match_Pose_Start_Count->m_iCount++;
		pNew = &poBase->m_pMatch_Pose[poB_Match_Pose_Start_Count->m_iStart + poB_Match_Pose_Start_Count->m_iCount];
	}
	pNew->m_iPose_B = iPose_A;
	pNew->m_oMatch_Index_Start_Count = oMatch_Index_Start_Count;

	return 1;
}

template int Pose_Base_Add_Feature(Pose_Base<double>* poBase, float Key_Point[][2], unsigned char Feature[][128], int iCount, int* piPose_ID);
template<typename _T>int Pose_Base_Add_Feature(Pose_Base<_T>* poBase, float Key_Point[][2], unsigned char Feature[][128],int iCount, int *piPose_ID)
{
	if (poBase->m_iKey_Point_Cur + iCount > poBase->m_iMax_Key_Point_Count)
	{
		printf("not inplemented yet\n");
		return 0;
	}
	memcpy(poBase->m_pKey_Point_Buffer + poBase->m_iKey_Point_Cur, Key_Point, iCount * 2 * sizeof(float));
	memcpy(poBase->m_pFeature_Buffer + poBase->m_iKey_Point_Cur, Feature, iCount * 128);
	poBase->m_pPose[poBase->m_iPose_Cur].m_oKey_Point_Start_Count = { (int)poBase->m_iKey_Point_Cur,(int)iCount };
	poBase->m_iKey_Point_Cur += iCount;
	
	*piPose_ID = poBase->m_iPose_Cur++;
	poBase->m_iPose_Count++;
	return 1;
}

void Transpose_8x8(
	unsigned char Source[], int iSource_Stride,
	unsigned char Dest[], int iDest_Stride)
{
	// 1. 读 8 行，每行 8 字节，放到 8 个 xmm 的低 64 位。
	__m128i r0 = _mm_loadl_epi64((const __m128i*)(Source + 0 * iSource_Stride));
	__m128i r1 = _mm_loadl_epi64((const __m128i*)(Source + 1 * iSource_Stride));
	__m128i r2 = _mm_loadl_epi64((const __m128i*)(Source + 2 * iSource_Stride));
	__m128i r3 = _mm_loadl_epi64((const __m128i*)(Source + 3 * iSource_Stride));
	__m128i r4 = _mm_loadl_epi64((const __m128i*)(Source + 4 * iSource_Stride));
	__m128i r5 = _mm_loadl_epi64((const __m128i*)(Source + 5 * iSource_Stride));
	__m128i r6 = _mm_loadl_epi64((const __m128i*)(Source + 6 * iSource_Stride));
	__m128i r7 = _mm_loadl_epi64((const __m128i*)(Source + 7 * iSource_Stride));

	// 2. 第 1 级：按字节交错相邻两行。
	__m128i t0 = _mm_unpacklo_epi8(r0, r1);
	__m128i t1 = _mm_unpacklo_epi8(r2, r3);
	__m128i t2 = _mm_unpacklo_epi8(r4, r5);
	__m128i t3 = _mm_unpacklo_epi8(r6, r7);

	// 3. 第 2 级：按 16 位交错。
	__m128i u0 = _mm_unpacklo_epi16(t0, t1);
	__m128i u1 = _mm_unpackhi_epi16(t0, t1);
	__m128i u2 = _mm_unpacklo_epi16(t2, t3);
	__m128i u3 = _mm_unpackhi_epi16(t2, t3);

	// 4. 第 3 级：按 32 位交错。
	__m128i v0 = _mm_unpacklo_epi32(u0, u2);
	__m128i v1 = _mm_unpackhi_epi32(u0, u2);
	__m128i v2 = _mm_unpacklo_epi32(u1, u3);
	__m128i v3 = _mm_unpackhi_epi32(u1, u3);

	// 5. 存储：每个 v 含两行，各 8 字节。
	_mm_storel_epi64((__m128i*)(Dest + 0 * iDest_Stride), v0);
	_mm_storel_epi64((__m128i*)(Dest + 1 * iDest_Stride), _mm_srli_si128(v0, 8));
	_mm_storel_epi64((__m128i*)(Dest + 2 * iDest_Stride), v1);
	_mm_storel_epi64((__m128i*)(Dest + 3 * iDest_Stride), _mm_srli_si128(v1, 8));
	_mm_storel_epi64((__m128i*)(Dest + 4 * iDest_Stride), v2);
	_mm_storel_epi64((__m128i*)(Dest + 5 * iDest_Stride), _mm_srli_si128(v2, 8));
	_mm_storel_epi64((__m128i*)(Dest + 6 * iDest_Stride), v3);
	_mm_storel_epi64((__m128i*)(Dest + 7 * iDest_Stride), _mm_srli_si128(v3, 8));
}

template<typename _T>void Reshape_in_Block_Col(_T Source[], int w, int h, int iBlock_Width, _T Dest[])
{
	/*for (int y = 0; y < h; y++)
		for (int x = 0; x < w; x++)
			Source[y * w + x] = x;*/

	const int iSub_Block_Width = 8;
	if (sizeof(_T) == 1)
	{//字节的玩法
		int y,h_In_Block = (h / iBlock_Width) * iBlock_Width;
		for (y = 0; y < h_In_Block; y += iBlock_Width)
		{
			for (int y1 = 0; y1 < iBlock_Width; y1 += iSub_Block_Width)
				for (int x = 0; x < w; x += iSub_Block_Width)
				{
					int iSource_Pos = (y + y1) * w + x;
					int iDest_Pos = y *128 +  x * iBlock_Width + y1;
					/*if (y == 32 && y1 == 0 && x == 0)
						printf("here");;*/

					//Transpose_32x32(&Source[iSource_Pos], w, &Dest[iDest_Pos], iBlock_Width);
					//Transpose_16x16(&Source[iSource_Pos], w, &Dest[iDest_Pos], iBlock_Width);
					Transpose_8x8(&Source[iSource_Pos], w, &Dest[iDest_Pos], iBlock_Width);
					/*if(y==32)
						Disp(&Dest[128*32], 32, 32, "Dest");*/
					//if (Dest[32] != 27)
						//printf("here");
				}
		}

		//余下照抄了事，保持横向
		int iPos = y * w;
		memcpy(&Dest[iPos], &Source[iPos], w * h - iPos);
	}
	//Disp(Source, 1, 8, "Source");
	//Disp(Dest, 32, 32, "Dest");
	return;
}
template float (*Pose_Base_Get_Key_Point(Pose_Base<double>* poPose, int iPose_ID, int* piCount))[2];
template<typename _T> float (*Pose_Base_Get_Key_Point(Pose_Base<_T>* poPose, int iPose_ID, int* piCount))[2]
{//根据Pose_ID 返回关键点
	Start_Count oStart_Count = poPose->m_pPose[iPose_ID].m_oKey_Point_Start_Count;
	if (piCount)
		*piCount = oStart_Count.m_iCount;
	return  &poPose->m_pKey_Point_Buffer[oStart_Count.m_iStart];
}

template unsigned char (*Pose_Base_Get_Feature(Pose_Base<double>* poPose, int iPose_ID, int* piCount))[128];
template<typename _T> unsigned char (*Pose_Base_Get_Feature(Pose_Base<_T>* poPose, int iPose_ID, int* piCount))[128]
{
	Start_Count oStart_Count = poPose->m_pPose[iPose_ID].m_oKey_Point_Start_Count;
	*piCount = oStart_Count.m_iCount;
	return  &poPose->m_pFeature_Buffer[oStart_Count.m_iStart];
}

template<typename _T>void Get_Feature_by_Block_Col(Pose_Base<_T> *poPose, int iPose_ID, unsigned char (**ppFeature),int *piCount)
{//目标以128行，32字节为一块
	Start_Count oStart_Count= poPose->m_pPose[iPose_ID].m_oKey_Point_Start_Count;
	*piCount = oStart_Count.m_iCount;
	
	unsigned char* pSource = poPose->m_pFeature_Buffer[oStart_Count.m_iStart];
	unsigned char* pFeature = (unsigned char*)pMalloc(128 * oStart_Count.m_iCount);
	if (!pFeature)
		*ppFeature = NULL;
	else
		*ppFeature = pFeature;

	unsigned long long tStart = iGet_Tick_Count();
	//for(int i=0;i<10000;i++)
	Reshape_in_Block_Col<unsigned char>(pSource, 128, oStart_Count.m_iCount, 32, pFeature);
	//printf("%lld\n", iGet_Tick_Count() - tStart);
	return;
}

static unsigned int iGet_Distance_AVX512(unsigned char A[], unsigned char  B[])
{
	//计算欧氏距离
#pragma pack(16)
	register __m512i Value_32;
	register __m512i Sum;
	union {
		__m512i Sum_32;
		__m256i Sum_16[2];
	};
#pragma pack()
	//此处还要改一改，改成128字节对齐，可能会快些
	Value_32 = _mm512_sub_epi16(_mm512_cvtepu8_epi16(*(__m256i*) & A[0]), _mm512_cvtepu8_epi16(*(__m256i*) & B[0]));
	Sum_32 = _mm512_mullo_epi16(Value_32, Value_32);

	Value_32 = _mm512_sub_epi16(_mm512_cvtepu8_epi16(*(__m256i*) & A[32]), _mm512_cvtepu8_epi16(*(__m256i*) & B[32]));
	Sum_32 = _mm512_adds_epu16(Sum_32, _mm512_mullo_epi16(Value_32, Value_32));
	Sum = _mm512_add_epi32(_mm512_cvtepu16_epi32(Sum_16[0]), _mm512_cvtepu16_epi32(Sum_16[1]));

	Value_32 = _mm512_sub_epi16(_mm512_cvtepu8_epi16(*(__m256i*) & A[64]), _mm512_cvtepu8_epi16(*(__m256i*) & B[64]));
	Sum_32 = _mm512_mullo_epi16(Value_32, Value_32);

	Value_32 = _mm512_sub_epi16(_mm512_cvtepu8_epi16(*(__m256i*) & A[96]), _mm512_cvtepu8_epi16(*(__m256i*) & B[96]));
	Sum_32 = _mm512_adds_epu16(Sum_32, _mm512_mullo_epi16(Value_32, Value_32));
	Sum = _mm512_add_epi32(Sum, _mm512_add_epi32(_mm512_cvtepu16_epi32(Sum_16[0]), _mm512_cvtepu16_epi32(Sum_16[1])));

	return _mm512_reduce_add_epi32(Sum);	//快些
}

#define Dot_1(A, B, iDistance) \
{ \
	unsigned char *A1=A,*B1=B;	\
	__m512i Sum;	\
	union {	\
		__m512i Sum_32;	\
		__m256i Sum_16[2];	\
	};	\
	Sum_32 = _mm512_mullo_epi16(_mm512_cvtepu8_epi16(*(__m256i*) & A1[0]), _mm512_cvtepu8_epi16(*(__m256i*) & B1[0]));	\
	Sum = _mm512_add_epi32(_mm512_cvtepu16_epi32(Sum_16[0]), _mm512_cvtepu16_epi32(Sum_16[1]));	\
	Sum_32 = _mm512_mullo_epi16(_mm512_cvtepu8_epi16(*(__m256i*) & A1[32]), _mm512_cvtepu8_epi16(*(__m256i*) & B1[32]));	\
	Sum = _mm512_add_epi32(Sum, _mm512_add_epi32(_mm512_cvtepu16_epi32(Sum_16[0]), _mm512_cvtepu16_epi32(Sum_16[1])));	\
	Sum_32 = _mm512_mullo_epi16(_mm512_cvtepu8_epi16(*(__m256i*) & A1[64]), _mm512_cvtepu8_epi16(*(__m256i*) & B1[64]));	\
	Sum = _mm512_add_epi32(Sum, _mm512_add_epi32(_mm512_cvtepu16_epi32(Sum_16[0]), _mm512_cvtepu16_epi32(Sum_16[1])));	\
	Sum_32 = _mm512_mullo_epi16(_mm512_cvtepu8_epi16(*(__m256i*) & A1[96]), _mm512_cvtepu8_epi16(*(__m256i*) & B1[96]));	\
	Sum = _mm512_add_epi32(Sum, _mm512_add_epi32(_mm512_cvtepu16_epi32(Sum_16[0]), _mm512_cvtepu16_epi32(Sum_16[1])));	\
	(iDistance)=_mm512_reduce_add_epi32(Sum);\
}

int iFind_Nearest(unsigned char A[128],int iIndex_A, unsigned char Feature[][128], int iCount, Neighbour_Item<unsigned int> Neighbour_B[][2])
{//完全无卵用，更慢
	int i;
	unsigned int iDistance;
	if (iCount < 2)
	{
		printf("Point Count less then 2\n");
		return -1;	//不干也罢
	}
	Neighbour_Item<unsigned int> Neighbour[2] = { -1,0xFFFFFFFF,-1,0xFFFFFFFF };
	for (i = 0; i < iCount; i++)
	{
		iDistance = iGet_Distance_AVX512(A, Feature[i]);
		if(iDistance<Neighbour[1].m_fDistance)
		{
			Neighbour[1] = { i,iDistance };
			if (Neighbour[1].m_fDistance < Neighbour[0].m_fDistance)
				swap(Neighbour[0], Neighbour[1]);
		}

		//Neighbour_B的第i点，也要判断一下是否为最近点
		if (iDistance < Neighbour_B[i][1].m_fDistance)
		{
			if (iDistance < Neighbour_B[i][0].m_fDistance)
			{
				Neighbour_B[i][1] = Neighbour_B[i][0];
				Neighbour_B[i][0] = { iIndex_A,iDistance };
			}
			else
				Neighbour_B[i][1] = { iIndex_A,iDistance };
		}
	}

	//第二部分，求夹角
	//将欧式距离改为点积，即比角度
	Dot_1(A, Feature[Neighbour[0].m_iIndex], Neighbour[0].m_fDistance);
	Dot_1(A, Feature[Neighbour[1].m_iIndex], Neighbour[1].m_fDistance);

	//再将距离大者放在[0],次之放在[1]
	//看cos(x)图像，夹角越小，函数值越大。点积越大，函数值越大
	if (Neighbour[0].m_fDistance < Neighbour[1].m_fDistance)
		swap(Neighbour[0], Neighbour[1]);

	const float kDistNorm = 1.0f / (512.0f * 512.0f);
	const float max_distance = 0.7f,
		max_ratio = 0.8f;

	float best_dist_normed = (float)acos(std::min(kDistNorm * Neighbour[0].m_fDistance, 1.0f));
	float second_best_dist_normed = (float)acos(std::min(kDistNorm * Neighbour[1].m_fDistance, 1.0f));

	//首先判断最优点是否大于阈值max_distance,如果大了就跳过，没毛病
	if (best_dist_normed > max_distance)
		return -1;
	//再判断最优点是否比次优点的0.8倍还大，表示最优点离此优点不能太远，
	//这个判断尚未明白。次优点不是来捣乱的吗？
	if (best_dist_normed >= max_ratio * second_best_dist_normed)
		return -1;


	//感觉此处还可以再进行逆向检查了事
	//将欧式距离改为点积
	Dot_1(Feature[Neighbour[0].m_iIndex], A, Neighbour[0].m_fDistance);
	Dot_1(Feature[Neighbour[1].m_iIndex], A, Neighbour[1].m_fDistance);

	return Neighbour[0].m_iIndex;
}

int iFind_Nearest(unsigned char A[128], unsigned char* pFeature_Block, unsigned char Feature[][128], int iCount, unsigned int Aux[])
{//完全无卵用，更慢
	if (iCount < 2)
	{
		printf("Invalid Count in iFind_Nearest\n");
		return -1;
	}

	unsigned int* pSum;
	if (Aux)
		pSum = Aux;
	else
		pSum = (unsigned int*)pMalloc(iCount * sizeof(unsigned int));
	int y,iCount_Align_32 = (iCount >> 5) << 5;

	unsigned char* pCur = pFeature_Block;
	for (y = 0; y < iCount_Align_32; y += 32)
	{//每次吐出32个Sum
		__m512i Sum_16[2] = {};
		for (int i = 0; i < 128; i++, pCur+=32)
		{
			union {
				__m512i c;
				__m256i c_16[2];
			};
			c = _mm512_sub_epi16(_mm512_set1_epi16(A[i]), _mm512_cvtepu8_epi16(*(__m256i*)pCur));
			c = _mm512_mullo_epi16(c, c);
			Sum_16[0] = _mm512_add_epi32(Sum_16[0], _mm512_cvtepu16_epi32(c_16[0]));
			Sum_16[1] = _mm512_add_epi32(Sum_16[1], _mm512_cvtepu16_epi32(c_16[1]));
			//printf("%d\n", Sum_16[0].m512i_i32[0]);
		}
		_mm512_storeu_epi32(&pSum[y], Sum_16[0]);
		_mm512_storeu_epi32(&pSum[y + 16], Sum_16[1]);
	}
	//余数
	for (; y < iCount; y++)
		pSum[y] = iGet_Distance_AVX512(A, &pFeature_Block[y * 128]);
	//Disp(pSum, 1, iCount, "Sum");
	 
	Neighbour_Item<unsigned int> Neighbour[2];
	if (pSum[0] <= pSum[1])
	{
		Neighbour[0] = { 0,pSum[0] };
		Neighbour[1] = { 1,pSum[1] };
	}
	else
	{
		Neighbour[0] = { 1,pSum[1] };
		Neighbour[1] = { 0,pSum[0] };
	}
	for (int i = 2; i < iCount; i++)
	{
		if (pSum[i] < Neighbour[1].m_fDistance)
		{
			Neighbour[1] = { i,pSum[i] };
			if (Neighbour[1].m_fDistance < Neighbour[0].m_fDistance)
				swap(Neighbour[0], Neighbour[1]);
		}
	}
	if (!Aux)
		Free(pSum);

	//将欧式距离改为点积，即比角度
	Dot_1(A, Feature[Neighbour[0].m_iIndex], Neighbour[0].m_fDistance);
	Dot_1(A, Feature[Neighbour[1].m_iIndex], Neighbour[1].m_fDistance);

	//再将距离大者放在[0],次之放在[1]
	//看cos(x)图像，夹角越小，函数值越大。点积越大，函数值越大
	if (Neighbour[0].m_fDistance < Neighbour[1].m_fDistance)
		swap(Neighbour[0], Neighbour[1]);


	const float kDistNorm = 1.0f / (512.0f * 512.0f);
	const float max_distance = 0.7f,
		max_ratio = 0.8f;

	//注意，此处是阀值
	float best_dist_normed = (float)acos(std::min(kDistNorm * Neighbour[0].m_fDistance, 1.0f));
	float second_best_dist_normed = (float)acos(std::min(kDistNorm * Neighbour[1].m_fDistance, 1.0f));

	//首先判断最优点是否大于阈值max_distance,如果大了就跳过，没毛病
	if (best_dist_normed > max_distance)		//限制了夹角
		return -1;
	//再判断最优点是否比次优点的0.8倍还大，表示最优点离此优点不能太远，
	//这个判断尚未明白。次优点不是来捣乱的吗？
	if (best_dist_normed >= max_ratio * second_best_dist_normed)
		return -1;

	return Neighbour[0].m_iIndex;
}

void _Sift_Match_2_Image(unsigned char Feature_A[][128], int iCount_A,
	unsigned char Feature_B[][128], int iCount_B,
	unsigned short (**ppMatch_Pair)[2], int* piCount)
{
	unsigned short (*pMatch_Pair)[2] = (unsigned short(*)[2])pMalloc(Max(iCount_A, iCount_B) * sizeof(unsigned short));
	Neighbour_Item<unsigned int>(*pNeighbour_B)[2] = (Neighbour_Item<unsigned int>(*)[2])pMalloc(iCount_B * 2 * sizeof(Neighbour_Item<unsigned int>));	//点集A的最近邻	
	memset(pNeighbour_B, 0xFF, iCount_B * 2 * sizeof(Neighbour_Item<unsigned int>));
	int iMatch_Count = 0;
	for (int i = 0; i < iCount_A; i++)
	{
		int iMatch_B = iFind_Nearest(Feature_A[i], i, Feature_B, iCount_B, pNeighbour_B);
		if (iMatch_B != -1)
		{
			pMatch_Pair[iMatch_Count][0] = i;
			pMatch_Pair[iMatch_Count][1] = iMatch_B;
			iMatch_Count++;
		}
	}

	//再找一次从B到A的查找
	int i, j;
	const float kDistNorm = 1.0f / (512.0f * 512.0f);
	const float max_distance = 0.7f,
		max_ratio = 0.8f;
	float best_dist_normed, second_best_dist_normed;

	for (i = j = 0; i < iMatch_Count; i++)
	{
		int iMatch_A = pMatch_Pair[i][0],
			iMatch_B = pMatch_Pair[i][1];
		Neighbour_Item<unsigned int>Neighbour_B[2] = { pNeighbour_B[iMatch_B][0],pNeighbour_B[iMatch_B][1] };

		//将欧式距离改为点积
		Dot_1(Feature_A[Neighbour_B[0].m_iIndex], Feature_B[iMatch_B], Neighbour_B[0].m_fDistance);
		Dot_1(Feature_A[Neighbour_B[1].m_iIndex], Feature_B[iMatch_B], Neighbour_B[1].m_fDistance);

		if (Neighbour_B[0].m_fDistance < Neighbour_B[1].m_fDistance)
			swap(Neighbour_B[0], Neighbour_B[1]);

		if (Neighbour_B[0].m_iIndex != iMatch_A)
			continue;

		best_dist_normed = (float)acos(std::min(kDistNorm * Neighbour_B[0].m_fDistance, 1.0f));
		//首先判断最优点是否大于阈值max_distance,如果大了就跳过，没毛病
		if (best_dist_normed > max_distance)
			continue;
		second_best_dist_normed = (float)acos(std::min(kDistNorm * Neighbour_B[1].m_fDistance, 1.0f));
		//再判断最优点是否比次优点的0.8倍还大，表示最优点离此优点不能太远，
		//这个判断尚未明白。次优点不是来捣乱的吗？
		if (best_dist_normed >= max_ratio * second_best_dist_normed)
			continue;

		pMatch_Pair[j][0] = iMatch_A;
		pMatch_Pair[j][1] = iMatch_B;
		j++;
	}

	*piCount = iMatch_Count = j;
	*ppMatch_Pair = pMatch_Pair;
	Free(pNeighbour_B);
	Shrink(pMatch_Pair, iMatch_Count * 2 * sizeof(unsigned short));
	return;
}



template<typename _T, int n>int Add_To_Hash_Table(Bin_Hash_Item<_T, n> Item[], int* piCur_Item, struct Bin_Hash_Item<_T, n> oItem, unsigned int Hash_Table[], int iHash_Size)
{
	unsigned int iPos, iCur_Item = *piCur_Item,
		iHash_Pos = iGet_Hash_Code(oItem.Bin, n) % iHash_Size;
	Bin_Hash_Item<_T, n> oItem_Exist;
	if ((iPos = Hash_Table[iHash_Pos]) != 0xFFFFFFFF)
	{
		do {
			oItem_Exist = Item[iPos];
			//if (oItem_Exist.m_iPos == oItem.m_iPos)
			if(memcmp(oItem_Exist.Bin,oItem.Bin,n)==0)//重复了，此点不加入散列表
				return 0;
		} while ((iPos = oItem_Exist.m_iNext) != 0xFFFFFFFF);
	}
	oItem.m_iNext = Hash_Table[iHash_Pos];
	Item[iCur_Item] = oItem;
	Hash_Table[iHash_Pos] = iCur_Item++;

//HASH_END:
//	;
	*piCur_Item = iCur_Item;
	return 1;
}

void Sift_Match_2_Image(float Key_Point_A[][2], unsigned char Feature_A[][128], int iCount_A,
	float Key_Point_B[][2], unsigned char Feature_B[][128], int iCount_B,
	unsigned short (**ppMatch_Pair)[2], int* piCount, float fPeak_Thresh)
{//这个才是有效函数，必须要把关键点考虑进来，去除重复点
	unsigned short(*pMatch_Pair)[2];
	int iMatch_Count;
	_Sift_Match_2_Image(Feature_A, iCount_A, Feature_B, iCount_B, &pMatch_Pair, &iMatch_Count);
	
	//Hash 查重
	const int n = sizeof(float) * 4;
	typedef Bin_Hash_Item<int, n> Hash_Item;

	unsigned int iHash_Size = iMatch_Count * 2 + 1;
	int iSize = ALIGN_SIZE_8(iHash_Size * sizeof(unsigned int)) +
		ALIGN_SIZE_8(iHash_Size * sizeof(Hash_Item));

	unsigned int* pHash_Table = (unsigned int*)pMalloc(iSize);
	memset(pHash_Table, -1, iHash_Size * sizeof(unsigned int));
	Hash_Item* pHash_Item = (Hash_Item*)ALIGN_SIZE_8(pHash_Table + iHash_Size);
	int i,j,iCur_Item = 0;
	for (i = j = 0; i < iMatch_Count; i++)
	{
		/*if (i >= 43)
			printf("here");*/
		Hash_Item oItem;
		memcpy(oItem.Bin, Key_Point_A[pMatch_Pair[i][0]], 2 * sizeof(float));
		memcpy(&oItem.Bin[2 * sizeof(float)], Key_Point_B[pMatch_Pair[i][1]], 2 * sizeof(float));
		oItem.m_iNext = -1;

		int iResult = Add_To_Hash_Table<int,n>(pHash_Item, &iCur_Item, oItem, pHash_Table, iHash_Size);
		if (iResult)
		{//成功加入
			pMatch_Pair[j][0] = pMatch_Pair[i][0];
			pMatch_Pair[j][1] = pMatch_Pair[i][1];
			j++;
		}else
		{
			//printf("Error:%d\n", i);
		}
	}
	if (piCount)
		*piCount = j;
	if (ppMatch_Pair)
		*ppMatch_Pair = pMatch_Pair;
	Free(pHash_Table);
	return;
}

template void Pose_Base_Get_Match_Key_Point(Pose_Base<double>* poBase, int iPose_A, int iPose_B, double(**ppPoint_A)[2], double(**ppPoint_B)[2], int* piCount);
template<typename _T1,typename _T2>void Pose_Base_Get_Match_Key_Point(Pose_Base<_T1>* poBase, int iPose_A, int iPose_B, _T2(**ppPoint_A)[2], _T2(**ppPoint_B)[2], int* piCount)
{//给定一对Pose_ID，返回两组匹配点
	Match_Pose* poMatch_Pose = Pose_Base_Get_Match_Pose<_T1>(poBase,iPose_A, iPose_B);
	int iMatch_Count = poMatch_Pose->m_oMatch_Index_Start_Count.m_iCount;
	if (!poMatch_Pose)
	{
		*ppPoint_A = *ppPoint_B = NULL;
		return;
	}
	_T2 (*pPoint_A)[2] = (_T2(*)[2])pMalloc(iMatch_Count * 2 * 2 * sizeof(_T2)),
		(* pPoint_B)[2] = pPoint_A + iMatch_Count;
	if (!pPoint_A)
	{
		*ppPoint_A = *ppPoint_B = NULL;
		return;
	}
	unsigned short (* pIndex)[2] = &poBase->m_pMatch_Index[poMatch_Pose->m_oMatch_Index_Start_Count.m_iStart];
	//Disp((unsigned short*)(pIndex + 43), 2, 2, "Index");

	float(* pA)[2] = &poBase->m_pKey_Point_Buffer[poBase->m_pPose[iPose_A].m_oKey_Point_Start_Count.m_iStart],
		(* pB)[2] = &poBase->m_pKey_Point_Buffer[poBase->m_pPose[iPose_B].m_oKey_Point_Start_Count.m_iStart];
	//Disp((float*)pA, poBase->m_pPose[iPose_A].m_oKey_Point_Start_Count.m_iCount, 2, "A");

	int iCount = poMatch_Pose->m_oMatch_Index_Start_Count.m_iCount;
	if(iPose_A<iPose_B)
	{
		for (int i = 0; i <iCount ; i++)
		{
			pPoint_A[i][0] = pA[pIndex[i][0]][0];
			pPoint_A[i][1] = pA[pIndex[i][0]][1];
			pPoint_B[i][0] = pB[pIndex[i][1]][0];
			pPoint_B[i][1] = pB[pIndex[i][1]][1];
		}
	}else
	{//要反转
		for (int i = 0; i < iCount; i++)
		{
			pPoint_A[i][0] = pA[pIndex[i][1]][0];
			pPoint_A[i][1] = pA[pIndex[i][1]][1];
			pPoint_B[i][0] = pB[pIndex[i][0]][0];
			pPoint_B[i][1] = pB[pIndex[i][0]][1];
		}
	}
	*ppPoint_A = pPoint_A;
	*ppPoint_B = pPoint_B;
	*piCount = iCount;
	//unsigned short (*pMatch_Index)[2] = poBase->m_pMatch_Index[]
}

static int Gauss_Filter_AVX512_1(float* pSource, int iWidth, int iHeight, int r, float* pDest,float *pAux_Buffer=NULL)
{//这个函数毫无意义，在写目标的时候相当于随机，所以分块还是分裂都更慢
	//记录一下避坑
	int iSize = ALIGN_SIZE_8(r * 2 + 1) +
		ALIGN_SIZE_16(iWidth * iHeight);

	float* pFilter = (float*)pMalloc(iSize * sizeof(float)),
		* pAux = pFilter + ALIGN_SIZE_8(r * 2 + 1);

	Gen_Gauss_Filter(r, pFilter);
	//Disp(pFilter, 1, r, "Filter");
	//Gen_Gauss_Filter(r, 0, pFilter);
	//memset(pFilter, 0, (r * 2 + 1) * sizeof(float));
	//pFilter[r] = 1;

	//第一步，纵向,全部是源尺寸
	short iBlock_Size_w = Min(iWidth, 16), iBlock_Size_h = Min(iHeight, 16),
		iW_Remain = iWidth & 0xF, iH_Remain = iHeight & 0xF;

	int iDest_Col_Size = iWidth * iBlock_Size_h;
	unsigned short	iBlock_Size_Sqr_Full = iBlock_Size_w * iBlock_Size_h,
		iBlock_Size_Sqr_Remain = iH_Remain * iBlock_Size_w,
		iBlock_Size_Sqr;

	float* pFilter_1 = &pFilter[r];
	__m512 Pixel_0_16, Pixel_1_16, Sum_16;
	__m512i Interval, Interval_Block, Interval_H_Remain/*, Interval_W_Remain*/;
	
	//目标为pAux，
	for (int i = 0; i < 16; i++)
	{
		Interval_Block.m512i_i32[i] = i * iBlock_Size_h;
		Interval_H_Remain.m512i_i32[i] = i * iH_Remain;;
	}
	__mmask16 Mask = (1 << iW_Remain) - 1;
	for (int x = 0; x < iWidth; x += iBlock_Size_w)
	{
		int iRemain_x = iWidth - x;
		if (iRemain_x >= iBlock_Size_w)
			iRemain_x = iBlock_Size_w;
		Interval = Interval_Block;
		iBlock_Size_Sqr = iBlock_Size_Sqr_Full;
		for (int y = 0; y < iHeight;)
		{
			int iRemain_y = iHeight - y;
			if (iRemain_y >= iBlock_Size_h)
				iRemain_y = iBlock_Size_h;
			else
			{
				Interval = Interval_H_Remain;
				iBlock_Size_Sqr = iBlock_Size_Sqr_Remain;
			}
			int iDest_Pos = (y >> 4) * iDest_Col_Size + (x >> 4) * iBlock_Size_Sqr;
			for (int iLine = 0; iLine < iRemain_y; iLine++, y++, iDest_Pos++)
			{
				//Sum_16 = _mm512_set1_ps(0.f);
				Pixel_0_16 = _mm512_loadu_ps(&pSource[x + iWidth * y]);
				Sum_16 = _mm512_add_ps(_mm512_set1_ps(0.f), _mm512_mul_ps(Pixel_0_16, _mm512_set1_ps(*pFilter_1)));

				for (int i = 1; i <= r; i++)
				{//以自己为中心，向上下两个方向逐步散出去
					int y1 = y - i;	//y1表示自己上面第i行，Pixel_0_16表示从内存中连续读入16个数据进寄存器
					Pixel_0_16 = _mm512_loadu_ps(&pSource[x + iWidth * (y1 < 0 ? 0 : y1)]);
					y1 = y + i; //y1表示下面第i行，Pixel_1_16表示从内存中连续读入16个数据进寄存器
					Pixel_1_16 = _mm512_loadu_ps(&pSource[x + iWidth * (y1 >= iHeight ? iHeight - 1 : y1)]);
					//Sum = (Pixel_0 + Pixel_1) * Filter[i]，此处利用了一个性质，以自己为中心，高斯核竖着看
					//上面第i个数据与限免第i个数据对应同一个权值，这是由高斯核的对成性决定的。
					Sum_16 = _mm512_add_ps(Sum_16, _mm512_mul_ps(_mm512_add_ps(Pixel_0_16, Pixel_1_16), _mm512_set1_ps(pFilter_1[i])));
				}
				//得到一样，马上可以转置，散到目标列上
				//if (iRemain_x == iBlock_Size_w)
				if (iRemain_x == 16)
					_mm512_i32scatter_ps(&pAux[iDest_Pos], Interval, Sum_16, 4);
				else
					_mm512_mask_i32scatter_ps(&pAux[iDest_Pos], Mask, Interval, Sum_16, 4);
			}
		}
	}

	//Disp(pAux, 17, 4, "B");
	//Disp(pAux + 17 * 16, 17, 2, "B");
	//Disp(pAux + 16 * 8*2, 2, 8, "B");

	//第二步，横向，由于已经转置，所以全部尺寸要改
	swap(iWidth, iHeight);
	swap(iBlock_Size_w, iBlock_Size_h);
	swap(iW_Remain, iH_Remain);
	iDest_Col_Size = iHeight * iBlock_Size_w;

	//目标为pAux，
	for (int i = 0; i < 16; i++)
	{//目标/写入 的块
		Interval_Block.m512i_i32[i] = i * iHeight;
		Interval_H_Remain.m512i_i32[i] = i * iH_Remain;
	}
	Mask = (1 << iW_Remain) - 1;

	float* pCur = pAux;
	for (int x = 0; x < iWidth; x += iBlock_Size_w)
	{//从左扫右，逐个Block 列扫过去
		int iRemain_x = iWidth - x;
		if (iRemain_x >= iBlock_Size_w)
			iRemain_x = iBlock_Size_w;
		Interval = Interval_Block;

		int iDest_Pos = x * iHeight;
		for (int y = 0; y < iHeight; y ++,iDest_Pos++)
		{
			Sum_16 = _mm512_mul_ps(_mm512_loadu_ps(&pCur[y * iRemain_x]), _mm512_set1_ps(*pFilter_1));
			for (int i = 1; i <= r; i++)
			{//以自己为中心，向上下两个方向逐步散出去
				int y1 = y - i;	//y1表示自己上面第i行，Pixel_0_16表示从内存中连续读入16个数据进寄存器
				Pixel_0_16 = _mm512_loadu_ps(&pCur[(y1 < 0 ? 0 : y1) * iRemain_x]);
				y1 = y + i; //y1表示下面第i行，Pixel_1_16表示从内存中连续读入16个数据进寄存器
				Pixel_1_16 = _mm512_loadu_ps(&pCur[(y1 >= iHeight ? iHeight - 1 : y1) * iRemain_x]);
				//上面第i个数据与限免第i个数据对应同一个权值，这是由高斯核的对成性决定的。
				Sum_16 = _mm512_add_ps(Sum_16, _mm512_mul_ps(_mm512_add_ps(Pixel_0_16, Pixel_1_16), _mm512_set1_ps(pFilter_1[i])));
			}
			if (iRemain_x == 16 /*&& iRemain_y == 16*/)
				_mm512_i32scatter_ps(&pDest[iDest_Pos], Interval, Sum_16, 4);
			else
				_mm512_mask_i32scatter_ps(&pDest[iDest_Pos], Mask, Interval, Sum_16, 4);
		}
		pCur += iRemain_x * iHeight;
	}

	//Disp(pDest, iWidth, iHeight, "B");
	Free(pFilter);
	return 1;
}

template int Sift_Match_2_Image(const char* pcFile_A, const char* pcFile_B, float(**ppPoint_A)[2], float(**ppPoint_B)[2], int* piCount, int o_min, float fPeak_Thresh);
template int Sift_Match_2_Image(const char* pcFile_A, const char* pcFile_B, double(**ppPoint_A)[2], double(**ppPoint_B)[2], int* piCount, int o_min, float fPeak_Thresh);
template<typename _T>int Sift_Match_2_Image(const char* pcFile_A, const char* pcFile_B, _T(**ppPoint_A)[2], _T(**ppPoint_B)[2], int* piCount, int o_min, float fPeak_Thresh)
{//两图匹配，找到Sift匹配掂，属于简易算法，绕过Pose Base

	float (*pKey_Point_A)[2] = NULL, (*pKey_Point_B)[2] = NULL;
	unsigned char (*pFeature_A)[128] = NULL, (*pFeature_B)[128] = NULL;
	unsigned short (*pMatch_Index)[2] = NULL;
	int iResult, iCount, iCount_A, iCount_B, bRet = 0;
	iResult = Get_Sift_Feature_1(pcFile_A, &pKey_Point_A, &pFeature_A, &iCount_A, o_min, fPeak_Thresh);
	iResult += Get_Sift_Feature_1(pcFile_B, &pKey_Point_B, &pFeature_B, &iCount_B, o_min, fPeak_Thresh);
	if (iResult < 2)goto END;

	Sift_Match_2_Image(pKey_Point_A, pFeature_A, iCount_A,
		pKey_Point_B, pFeature_B, iCount_B, &pMatch_Index, &iCount);
	if (!iCount)goto END;

	_T(*pPoint_A)[2], (*pPoint_B)[2];
	pPoint_A = (_T(*)[2])pMalloc(iCount * 2 * 2 * sizeof(_T));
	pPoint_B = pPoint_A + iCount;
	if (!pPoint_A)goto END;
	for (int i = 0; i < iCount; i++)
	{
		pPoint_A[i][0] = pKey_Point_A[pMatch_Index[i][0]][0];
		pPoint_A[i][1] = pKey_Point_A[pMatch_Index[i][0]][1];
		pPoint_B[i][0] = pKey_Point_B[pMatch_Index[i][1]][0];
		pPoint_B[i][1] = pKey_Point_B[pMatch_Index[i][1]][1];
		//printf("i:%d %f %f %f %f\n", i, pPoint_A[i][0], pPoint_A[i][1], pPoint_B[i][0], pPoint_B[i][1]);
	}
	*ppPoint_A = pPoint_A;
	*ppPoint_B = pPoint_B;
	if (piCount)
		*piCount = iCount;
	bRet = 1;
END:
	Free(pKey_Point_A);
	Free(pKey_Point_B);
	Free(pFeature_A);
	Free(pFeature_B);
	Free(pMatch_Index);
	return bRet;
}