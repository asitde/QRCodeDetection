#include "stdafx.h"
#include "QRCodeDetector.h"


const double QRCodeDetector::PIXEL_TO_MM = 0.09;

QRCodeDetector::QRCodeDetector()
	:m_strImagePath(""),
	m_strRootDirectory(""),
	m_dParam1(120),
	m_dParam2(60),
	m_iMinRadius(650), m_iMaxRadius(765)
{}

QRCodeDetector::~QRCodeDetector() 
{
	if (m_ofs.is_open()) {
		m_ofs.close();
	}
}

void QRCodeDetector::ApplyGaussianBlur(const cv::Mat& rOriginalImage, cv::Mat& rResultantImage)
{
	GaussianBlur(rOriginalImage, rResultantImage, cv::Size(9, 9), 2, 2); // Alternatively use 9,9
}

void QRCodeDetector::ApplyMedianBlur(const cv::Mat& rOriginalImage, cv::Mat& rResultantImage)
{
	medianBlur(rOriginalImage, rResultantImage, 5);
}

void QRCodeDetector::ConvertToGrayScale(const cv::Mat& rOriginalImage, cv::Mat& rResultantImage)
{
	cvtColor(rOriginalImage, rResultantImage, CV_BGR2GRAY);
}

void QRCodeDetector::ApplyThresold(const cv::Mat& rOriginalImage, cv::Mat& rResultantImage)
{                      
	threshold(rOriginalImage, rResultantImage, 128, 255, CV_THRESH_BINARY_INV);
}

void QRCodeDetector::ApplyAdaptiveThresold(const cv::Mat& rOriginalImage, cv::Mat& rResultantImage)
{
	adaptiveThreshold(rOriginalImage, rResultantImage, 255, cv::ADAPTIVE_THRESH_MEAN_C, cv::THRESH_BINARY, 3, 0);
}

void QRCodeDetector::ShowImage(const cv::Mat& rSrcImage, const std::string& rstrWindowName)
{
	////cv::namedWindow(rstrWindowName, CV_WINDOW_NORMAL);
	////cv::imshow(rstrWindowName, rSrcImage);
	////cv::waitKey(0);
	////cv::destroyWindow(rstrWindowName);
}



void QRCodeDetector::DrawDetectedCircles(cv::Mat& rOriginalImage, const std::vector<cv::Vec3f>& rVecCircles)
{
	////for (size_t i = 0; i < rVecCircles.size(); i++) {
	////	cv::Vec3i c = rVecCircles[i];
	////	circle(rOriginalImage, cv::Point(c[0], c[1]), c[2], cv::Scalar(0, 0, 255), 3, cv::LINE_AA);
	////	circle(rOriginalImage, cv::Point(c[0], c[1]), 2, cv::Scalar(255, 0, 0), 3, cv::LINE_AA);
	////}

	for (size_t i = 0; i < rVecCircles.size(); i++) {
		cv::Point center(cvRound(rVecCircles[i][0]), cvRound(rVecCircles[i][1]));
		int radius = cvRound(rVecCircles[i][2]);
		// circle center
		circle(rOriginalImage, center, 3, cv::Scalar(0, 255, 0), -1, 8, 0);
		// circle outline
		circle(rOriginalImage, center, radius, cv::Scalar(0, 0, 255), 3, 8, 0);
	}
}

void QRCodeDetector::SearchContours(const cv::Mat& rSrcImage, std::vector<std::vector<cv::Point>>& rVecContours)
{

	std::vector<std::vector<cv::Point> > contours;
	std::vector<cv::Vec4i> hierarchy;


	////cv::Mat canny_output;
	////cv::RNG rng;

	/////// Detect edges using canny
	////Canny(rSrcImage, canny_output, 1, 3, 3);
	/////// Find contours
	////findContours(canny_output, rVecContours, hierarchy, CV_RETR_TREE, CV_CHAIN_APPROX_SIMPLE, cv::Point(0, 0));

	/////// Draw contours
	////cv::Mat drawing = cv::Mat::zeros(canny_output.size(), CV_8UC3);
	////for (int i = 0; i< rVecContours.size(); i++)
	////{
	////	cv::Scalar color = cv::Scalar(rng.uniform(0, 255), rng.uniform(0, 255), rng.uniform(0, 255));
	////	drawContours(drawing, rVecContours, i, color, 2, 8, hierarchy, 0, cv::Point());
	////}

	////ShowImage(drawing, "Contours");

	findContours(rSrcImage, rVecContours, hierarchy, CV_RETR_LIST, CV_CHAIN_APPROX_SIMPLE, cv::Point(0, 0));



}


namespace nsUtility {

	static void RecursiveDelete(CString szPath)
	{
		CFileFind ff;
		CString path = szPath;

		if (path.Right(1) != "\\")
			path += "\\";

		path += "*.*";

		BOOL res = ff.FindFile(path);

		while (res)
		{
			res = ff.FindNextFile();
			if (!ff.IsDots() && !ff.IsDirectory())
				DeleteFile(ff.GetFilePath());
			else if (ff.IsDirectory())
			{
				path = ff.GetFilePath();
				////RecursiveDelete(path);
				RemoveDirectory(path);
			}
		}
		::RemoveDirectory(szPath);

	}
}

void QRCodeDetector::SetParams(double dParam1, double dParam2)
{
	m_dParam1 = dParam1;
	m_dParam2 = dParam2;
}

void QRCodeDetector::SetRadii(int iMinRadius, int iMaxRadius)
{
	m_iMinRadius = iMinRadius;
	m_iMaxRadius = iMaxRadius;
}

bool QRCodeDetector::Run(const CString& cstrImageDirectory)
{
	bool	bRetVal{ false };
	HANDLE	hFile{ NULL };


	while (true) {

		CString	cstrSearchDirectory{ cstrImageDirectory + _T("\\") };
		CString cstrRoot{ cstrSearchDirectory + _T("*")};

		WIN32_FIND_DATA file;
		hFile = FindFirstFile(cstrRoot, &file);
		if (hFile == INVALID_HANDLE_VALUE) {
			break;
		}


		// House keeping
		std::string strRootDirectory{ CT2CA{ cstrSearchDirectory }.m_psz };
		m_strRootDirectory = strRootDirectory;
		m_strOutputDirectory = m_strRootDirectory + "Output\\";

		// Create output directory
		CString cstrOutputDirectory{ cstrSearchDirectory + _T("Output") };
		nsUtility::RecursiveDelete(cstrOutputDirectory);
		::Sleep(3000);
		::CreateDirectory(cstrOutputDirectory,nullptr);
		std::string strCSV{ m_strOutputDirectory + "Output.csv"};
		m_ofs.open(strCSV, std::ofstream::out);
		m_ofs << "Image File" << "," << "Left" << "," << "Top" << "," << "Right" << "," << "Bottom" << std::endl;


		while (FindNextFile(hFile, &file)) {

			if (file.dwFileAttributes == FILE_ATTRIBUTE_DIRECTORY) {
				continue;
			}

			CString cstrImageFileName = file.cFileName;
			std::string strImageFileName{ CT2CA{ cstrImageFileName}.m_psz };
			//Process(strImageFileName);
			ProcessEx(strImageFileName);
		}

		bRetVal = true;

		break;
	}

	if ( (hFile != NULL) || (hFile != INVALID_HANDLE_VALUE) ) {
		::FindClose(hFile);
		hFile = NULL;
	}

	return bRetVal;
}

void QRCodeDetector::DetectCircles(const cv::Mat& rSrcImage, std::vector<cv::Vec3f>& rVecCircles)
{
	//HoughCircles(rSrcImage, rVecCircles, cv::HOUGH_GRADIENT, 1, rSrcImage.rows / 2, 100, 100, rSrcImage.rows/ 6, rSrcImage.rows/2);
	//HoughCircles(rSrcImage, rVecCircles, cv::HOUGH_GRADIENT, 1, rSrcImage.rows / 2, 100, 100, 0, 900);
	//HoughCircles(rSrcImage, rVecCircles, cv::HOUGH_GRADIENT, 1, rSrcImage.rows / 2, 200, 100, 0, 800);
	//HoughCircles(rSrcImage, rVecCircles, cv::HOUGH_GRADIENT, 1, rSrcImage.rows / 2, 200, 50, rSrcImage.rows/4, rSrcImage.rows);


	//HoughCircles(rSrcImage, rVecCircles, cv::HOUGH_GRADIENT, 1, rSrcImage.rows, 100, 60, rSrcImage.cols / 4, rSrcImage.rows / 4.5);
	////HoughCircles(rSrcImage, rVecCircles, cv::HOUGH_GRADIENT, 1, rSrcImage.rows, 60, 60, 650, 900);
	//HoughCircles(rSrcImage, rVecCircles, cv::HOUGH_GRADIENT, 1, rSrcImage.rows, 120, 60, 650, 775);
	////HoughCircles(rSrcImage, rVecCircles, cv::HOUGH_GRADIENT, 1, rSrcImage.rows, 120, 60, 650, 765);

	HoughCircles(rSrcImage, rVecCircles, cv::HOUGH_GRADIENT, 1, rSrcImage.rows, m_dParam1, m_dParam2, m_iMinRadius, m_iMaxRadius);


}


bool QRCodeDetector::ProcessEx(const std::string& rstrImageFile)
{
	bool bRetVal{ false };

	IplImage*		imgGrayScale = nullptr;
	CvMemStorage*	storage = nullptr;
	double			dRadius = 0;

	while (true) {

		std::string	strImagePath{ m_strRootDirectory + rstrImageFile };
		std::string	strOutputImagePath{ m_strOutputDirectory + rstrImageFile };


		cv::Mat originalImage = cv::imread(strImagePath);
		if (originalImage.empty()) {
			break;
		}

		//// No need to read the image as grayscale. We will read the imag as it is 
		//// and convert it to graysale after the necessary noise reduction.
		//cv::Mat grayScaleImage = cv::imread(strImagePath, CV_LOAD_IMAGE_GRAYSCALE);
		//if (grayScaleImage.empty()) {
		//	break;
		//}

		cv::Mat processedImage;
		ApplyMedianBlur(originalImage, processedImage);
		ShowImage(processedImage, strImagePath);
		ConvertToGrayScale(processedImage, processedImage);
		ShowImage(processedImage, strImagePath);

		std::vector<cv::Vec3f> vecCircles;
		DetectCircles(processedImage, vecCircles);
		if (vecCircles.size() == 0) {
			imwrite(strOutputImagePath, originalImage);
			m_ofs << rstrImageFile << "," << 0 << "," << 0 << "," << 0 << "," << 0 << "," << 0 << std::endl;
			break;
		}
		DrawDetectedCircles(originalImage, vecCircles);
		ShowImage(originalImage, strImagePath);


		dRadius = vecCircles[0][2];
		cv::Point ptCenter;
		ptCenter.x = vecCircles[0][0];
		ptCenter.y = vecCircles[0][1];

		int iRows = originalImage.rows;
		int iCols = originalImage.cols;

		imgGrayScale = cvLoadImage(strImagePath.c_str(), CV_LOAD_IMAGE_GRAYSCALE);
		//cvCanny(imgGrayScale, imgGrayScale, 0, thresh, 5);
		
		//cvCanny(imgGrayScale, imgGrayScale, 50, 200, 3);
		if (iRows > 1000) {
			cvSmooth(imgGrayScale, imgGrayScale, CV_MEDIAN, 9, 9, 0);
		}
		//cvDilate(imgGrayScale, imgGrayScale);
		//thresholding the grayscale image to get better results
		//cvThreshold(imgGrayScale, imgGrayScale, 128, 255, CV_THRESH_BINARY);

		cvAdaptiveThreshold(imgGrayScale, imgGrayScale, 255, CV_ADAPTIVE_THRESH_MEAN_C, CV_THRESH_BINARY , 13, 1);
		//converting the original image into grayscale
		//IplImage imgGrayScale(im); // cvCreateImage(cvGetSize(&crop), 8, 1);

		//thresholding the grayscale image to get better results
		//cvThreshold(&imgGrayScale, &imgGrayScale, 128, 255, CV_THRESH_BINARY);
		CvSeq*			contours = nullptr;  //hold the pointer to a contour in the memory block
		CvSeq*			result = nullptr;   //hold sequence of points of a contour
		storage = cvCreateMemStorage(0); //storage area for all contours

														 //finding all contours in the image
		cvFindContours(imgGrayScale, storage, &contours, sizeof(CvContour), CV_RETR_LIST, CV_CHAIN_APPROX_SIMPLE, cvPoint(0, 0));
		cv::Point ptCtrBoxQR[10];
		double dArmOR[10];
		int iQR = 0;
		int iYDiff = 0;
		int iXDiff = 0;
		int iTopDistance = 0;
		int iLeftDistance = 0;
		int iRightDistance = 0;
		int iBottomDistance = 0;
		double dArm = 0;
		double dArmMax = 0;
		double dArmMin = 9999;

		//iterating through each contour
		while (contours)
		{			//obtain a sequence of points of contour, pointed by the variable 'contour'
			result = cvApproxPoly(contours, sizeof(CvContour), storage, CV_POLY_APPROX_DP, cvContourPerimeter(contours)*0.02, 0);

			//if there are 4 vertices in the contour(It should be a quadrilateral)
			CvPoint *pt[4];
			cv::Point ptCtrBox;
			ptCtrBox.x = -1;
			ptCtrBox.y = -1;

			CvSeq* approx;
			approx = cvApproxPoly(contours, sizeof(CvContour), storage, CV_POLY_APPROX_DP, cvContourPerimeter(contours) * 0.02, 0);
		//	approx = cvApproxPoly(contours, sizeof(CvContour), storage, CV_POLY_APPROX_DP, cvarcLength(Mat(contours[i]), true)*0.02, 1);

				double d = abs(cvContourArea(approx, CV_WHOLE_SEQ, 0));

				int iMinArea = iRows / 4; // 100/ 2000
				int iMaxArea = iMinArea * 10; // 25000

			if((approx->total == 4 && abs(cvContourArea(approx, CV_WHOLE_SEQ, 0)) > iMinArea && abs(cvContourArea(approx, CV_WHOLE_SEQ, 0)) < iMaxArea && cvCheckContourConvexity(approx) != 0) && (result->total == 4))
			{
				//iterating through each point
				for (int i = 0; i < 4; i++) {
					pt[i] = (CvPoint*)cvGetSeqElem(result, i);
				}

				double dPtDistance = sqrt(
					(((*pt[0]).x - ptCenter.x)* ((*pt[0]).x - ptCenter.x)) +
					(((*pt[0]).y - ptCenter.y)* ((*pt[0]).y - ptCenter.y))
				);
				if (dPtDistance < dRadius) {
					if ((((abs((*pt[0]).x - (*pt[1]).x) > 15) ||
						(abs((*pt[0]).y - (*pt[1]).y) > 15)) &&
						((abs((*pt[1]).x - (*pt[2]).x) > 15) ||
						(abs((*pt[1]).y - (*pt[2]).y) > 15)) &&
							((abs((*pt[2]).x - (*pt[3]).x) > 15) ||
						(abs((*pt[2]).y - (*pt[3]).y) > 15)) &&
								((abs((*pt[3]).x - (*pt[0]).x) > 15) ||
						(abs((*pt[3]).y - (*pt[0]).y) > 15))) &&

									(((abs((*pt[0]).x - (*pt[1]).x) < 50) ||
						(abs((*pt[0]).y - (*pt[1]).y) < 50)) &&
										((abs((*pt[1]).x - (*pt[2]).x) < 50) ||
										(abs((*pt[1]).y - (*pt[2]).y) < 50)) &&
											((abs((*pt[2]).x - (*pt[3]).x) < 50) ||
										(abs((*pt[2]).y - (*pt[3]).y) < 50)) &&
												((abs((*pt[3]).x - (*pt[0]).x) < 50) ||
										(abs((*pt[3]).y - (*pt[0]).y) < 50))))
					{

						double dArm1 = sqrt(
							(((*pt[0]).x - (*pt[1]).x)* ((*pt[0]).x - (*pt[1]).x)) +
							(((*pt[0]).y - (*pt[1]).y)* ((*pt[0]).y - (*pt[1]).y))
						);
						double dArm2 = sqrt(
							(((*pt[1]).x - (*pt[2]).x)* ((*pt[1]).x - (*pt[2]).x)) +
							(((*pt[1]).y - (*pt[2]).y)* ((*pt[1]).y - (*pt[2]).y))
						);
						double dArm3 = sqrt(
							(((*pt[2]).x - (*pt[3]).x)* ((*pt[2]).x - (*pt[3]).x)) +
							(((*pt[2]).y - (*pt[3]).y)* ((*pt[2]).y - (*pt[3]).y))
						);
						double dArm4 = sqrt(
							(((*pt[3]).x - (*pt[0]).x)* ((*pt[3]).x - (*pt[0]).x)) +
							(((*pt[3]).y - (*pt[0]).y)* ((*pt[3]).y - (*pt[0]).y))
						);
						dArm = (dArm1 + dArm2 + dArm3 + dArm4) / 4;
						if (dArmMax < dArm)
							dArmMax = dArm;
						if (dArmMin > dArm)
							dArmMin = dArm;
						line(originalImage, *pt[0], *pt[1], cvScalar(0, 255, 0), 4);
						line(originalImage, *pt[1], *pt[2], cvScalar(0, 255, 0), 4);
						line(originalImage, *pt[2], *pt[3], cvScalar(0, 255, 0), 4);
						line(originalImage, *pt[3], *pt[0], cvScalar(0, 255, 0), 4);

						ptCtrBox.x = double((*pt[0]).x + (*pt[1]).x + (*pt[2]).x + (*pt[3]).x) / 4;
						ptCtrBox.y = double((*pt[0]).y + (*pt[1]).y + (*pt[2]).y + (*pt[3]).y) / 4;

						iYDiff = abs(ptCtrBox.y - (*pt[0]).y);
						iXDiff = abs(ptCtrBox.x - (*pt[0]).x);

						bool bExist = false;
						if (iQR > 0) {
							for (int iQ = 0; iQ < iQR; iQ++) {
								if ((abs(ptCtrBoxQR[iQ].x - ptCtrBox.x) < 5) && (abs(ptCtrBoxQR[iQ].y - ptCtrBox.y) < 5)) {
									bExist = true;
								}
							}
						}

						if (!bExist) {
							ptCtrBoxQR[iQR].x = ptCtrBox.x;
							ptCtrBoxQR[iQR].y = ptCtrBox.y;
							dArmOR[iQR] = dArm;
							iQR = iQR + 1;
						}
					}
				}
			}


			if ((dRadius != -1) && (ptCtrBox.x != -1) && (ptCtrBox.y != -1)) {
				double dDistanceBoxFromCtr = sqrt(((ptCtrBox.x - ptCenter.x) * (ptCtrBox.x - ptCenter.x)) +
					((ptCtrBox.y - ptCenter.y) * (ptCtrBox.y - ptCenter.y)));

				double dDistanceFromOutside = dRadius - dDistanceBoxFromCtr;

				line(originalImage, ptCtrBox, ptCenter, cvScalar(0, 255, 0), 2);

				//	cv::waitKey(0);
				int i = 0;
				i++;
			}

			contours = contours->h_next;
		}

		if (iQR >= 3) {
			if (iQR > 3) { // elemenate south west square
				while (true) {
					int iS1 = -1;
					int iS2 = -1;
					double dMax1 = 0;

					//if((ptCtrBoxQR[0].y > ptCtrBoxQR[1].y)
					for (int iE = 0; iE < iQR; iE++) {
						if (dMax1 < ptCtrBoxQR[iE].y) {
							dMax1 = ptCtrBoxQR[iE].y;
							iS1 = iE;
						}
					}

					dMax1 = 0;

					//if((ptCtrBoxQR[0].y > ptCtrBoxQR[1].y)
					for (int iE = 0; iE < iQR; iE++) {
						if (iE != iS1) {
							if (dMax1 < ptCtrBoxQR[iE].y) {
								dMax1 = ptCtrBoxQR[iE].y;
								iS2 = iE;
							}
						}
					}

					int iSW = -1;
					if ((iS1 >= 0) && (iS2 >= 0)) {
						if (ptCtrBoxQR[iS1].x > ptCtrBoxQR[iS2].x) {
							iSW = iS1;
						}
						else {
							iSW = iS2;
						}
					}

					if (iSW >= 0) {
						ptCtrBoxQR[iSW].x = -1;
						ptCtrBoxQR[iSW].y = -1;
						for (int iE = iSW; iE < iQR; iE++) {
							ptCtrBoxQR[iE].x = ptCtrBoxQR[iE + 1].x;
							ptCtrBoxQR[iE].y = ptCtrBoxQR[iE + 1].y;
							if (iE + 1 < 10) {
								ptCtrBoxQR[iE + 1].x = -1;
								ptCtrBoxQR[iE + 1].y = -1;
							}
						}
						iQR--;
					}
					if (iQR == 3) {
						break;
					}
				}
			}


			int iBottom = -1;
			int iLeft = -1;
			int iRight = -1;
			int iTop = -1;

			cv::Point ptLeft1;
			cv::Point ptTop1;
			cv::Point ptRight1;
			cv::Point ptBottom1;

			ptLeft1.x = -1;
			ptLeft1.y = -1;

			ptTop1.x = -1;
			ptTop1.y = -1;

			ptRight1.x = -1;
			ptRight1.y = -1;

			ptBottom1.x = -1;
			ptBottom1.y = -1;

			if ((ptCtrBoxQR[2].x < ptCtrBoxQR[1].x) && (ptCtrBoxQR[2].x < ptCtrBoxQR[0].x)) {
				iLeft = 2;
			}
			else if ((ptCtrBoxQR[1].x < ptCtrBoxQR[0].x) && (ptCtrBoxQR[1].x < ptCtrBoxQR[2].x)) {
				iLeft = 1;
			}
			else {
				iLeft = 0;
			}

			if ((ptCtrBoxQR[2].x > ptCtrBoxQR[1].x) && (ptCtrBoxQR[2].x > ptCtrBoxQR[0].x)) {
				iRight = 2;
			}
			else if ((ptCtrBoxQR[1].x > ptCtrBoxQR[0].x) && (ptCtrBoxQR[1].x > ptCtrBoxQR[2].x)) {
				iRight = 1;
			}
			else {
				iRight = 0;
			}

			if ((ptCtrBoxQR[2].y < ptCtrBoxQR[1].y) && (ptCtrBoxQR[2].y < ptCtrBoxQR[0].y)) {
				iTop = 2;
			}
			else if ((ptCtrBoxQR[1].y < ptCtrBoxQR[0].y) && (ptCtrBoxQR[1].y < ptCtrBoxQR[2].y)) {
				iTop = 1;
			}
			else {
				iTop = 0;
			}

			if ((ptCtrBoxQR[2].y > ptCtrBoxQR[1].y) && (ptCtrBoxQR[2].y > ptCtrBoxQR[0].y)) {
				iBottom = 2;
			}
			else if ((ptCtrBoxQR[1].y > ptCtrBoxQR[0].y) && (ptCtrBoxQR[1].y > ptCtrBoxQR[2].y)) {
				iBottom = 1;
			}
			else {
				iBottom = 0;
			}
			///*if ((abs(ptCtrBoxQR[2].y - ptCtrBoxQR[0].y) > dArm * 2) && (abs(ptCtrBoxQR[2].y - ptCtrBoxQR[1].y) > dArm * 2) &&
			//	(abs(ptCtrBoxQR[1].y - ptCtrBoxQR[0].y) < dArm * 2)) {
			//	iBottom = 2;

			//	if (ptCtrBoxQR[1].x > ptCtrBoxQR[0].x) {
			//		iLeft = 0;
			//		iRight = 1;
			//	}
			//	else {
			//		iLeft = 1;
			//		iRight = 0;
			//	}
			//	if (ptCtrBoxQR[1].y > ptCtrBoxQR[0].y) {
			//		iTop = 0;
			//	}
			//	else {
			//		iTop = 1;
			//	}

			//}
			//else if ((abs(ptCtrBoxQR[1].y - ptCtrBoxQR[0].y) > dArm * 2) && (abs(ptCtrBoxQR[1].y - ptCtrBoxQR[1].y) > dArm * 2) &&
			//	(abs(ptCtrBoxQR[2].y - ptCtrBoxQR[0].y) < dArm * 2)) {
			//	iBottom = 1;
			//	if (ptCtrBoxQR[2].x > ptCtrBoxQR[0].x) {
			//		iLeft = 0;
			//		iRight = 2;
			//	}
			//	else {
			//		iLeft = 2;
			//		iRight = 0;
			//	}
			//	if (ptCtrBoxQR[2].y > ptCtrBoxQR[0].y) {
			//		iTop = 0;
			//	}
			//	else {
			//		iTop = 1;
			//	}

			//}
			//else if ((abs(ptCtrBoxQR[0].y - ptCtrBoxQR[1].y) > dArm * 2) && (abs(ptCtrBoxQR[0].y - ptCtrBoxQR[2].y) > dArm * 2) &&
			//	(abs(ptCtrBoxQR[2].y - ptCtrBoxQR[1].y) < dArm * 2)) {
			//	iBottom = 0;
			//	if (ptCtrBoxQR[2].x > ptCtrBoxQR[1].x) {
			//		iLeft = 1;
			//		iRight = 2;
			//	}
			//	else {
			//		iLeft = 2;
			//		iRight = 1;
			//	}
			//	if (ptCtrBoxQR[1].y > ptCtrBoxQR[1].y) {
			//		iTop = 1;
			//	}
			//	else {
			//		iTop = 2;
			//	}

			//}*/
			if (iLeft > -1 && iRight > -1 && iBottom > -1 && iTop > -1) {
				if ((abs(ptCtrBoxQR[iLeft].y - ptCtrBoxQR[iRight].y)) < 25) {
					ptTop1.x = (ptCtrBoxQR[iLeft].x + ptCtrBoxQR[iRight].x) / 2;
					ptTop1.y = (ptCtrBoxQR[iLeft].y + ptCtrBoxQR[iRight].y) / 2;
				}
				else if (ptCtrBoxQR[iLeft].y < ptCtrBoxQR[iRight].y) {
					ptTop1.x = ptCtrBoxQR[iLeft].x;
					ptTop1.y = ptCtrBoxQR[iLeft].y;
				}
				else {
					ptTop1.x = ptCtrBoxQR[iRight].x;
					ptTop1.y = ptCtrBoxQR[iRight].y;
				}

				if ((abs(ptCtrBoxQR[iLeft].x - ptCtrBoxQR[iBottom].x)) < 25) {
					ptLeft1.x = (ptCtrBoxQR[iLeft].x + ptCtrBoxQR[iBottom].x) / 2;
					ptLeft1.y = (ptCtrBoxQR[iLeft].y + ptCtrBoxQR[iBottom].y) / 2;
				}
				else if (ptCtrBoxQR[iLeft].x < ptCtrBoxQR[iBottom].x) {
					ptLeft1.x = ptCtrBoxQR[iLeft].x;
					ptLeft1.y = ptCtrBoxQR[iLeft].y;
				}

				if ((abs(ptCtrBoxQR[iRight].x - ptCtrBoxQR[iBottom].x)) < 25) {
					ptRight1.x = (ptCtrBoxQR[iRight].x + ptCtrBoxQR[iBottom].x) / 2;
					ptRight1.y = (ptCtrBoxQR[iRight].y + ptCtrBoxQR[iBottom].y) / 2;
				}
				else if (ptCtrBoxQR[iRight].x > ptCtrBoxQR[iBottom].x) {
					ptRight1.x = ptCtrBoxQR[iRight].x;
					ptRight1.y = ptCtrBoxQR[iRight].y;
				}
				/*else {
				ptRight.x = ptCtrBoxQR[2].x;
				ptRight.y = ptCtrBoxQR[2].y;
				}*/

				if ((abs(ptCtrBoxQR[iBottom].y - ptCtrBoxQR[iRight].y)) < 25) {
					ptBottom1.x = (ptCtrBoxQR[iBottom].x + ptCtrBoxQR[iRight].x) / 2;
					ptBottom1.y = (ptCtrBoxQR[iBottom].y + ptCtrBoxQR[iRight].y) / 2;
				}
				else if ((abs(ptCtrBoxQR[iBottom].y - ptCtrBoxQR[iLeft].y)) < 25) {
					ptBottom1.x = (ptCtrBoxQR[iBottom].x + ptCtrBoxQR[iLeft].x) / 2;
					ptBottom1.y = (ptCtrBoxQR[iBottom].y + ptCtrBoxQR[iLeft].y) / 2;
				}
				else {
					ptBottom1.x = ptCtrBoxQR[iBottom].x;
					ptBottom1.y = ptCtrBoxQR[iBottom].y;
				}
			}

			cv::Point ptLeft;
			cv::Point ptRight;
			cv::Point ptBottom;
			cv::Point ptTop;

			//ptLeft.y = ptLeft1.y;
			//ptRight.y = ptRight1.y;
			//ptBottom.x = ptBottom1.x;
			//ptTop.x = ptTop1.x;

			ptLeft.y = ptCtrBoxQR[iLeft].y;
			ptRight.y = ptCtrBoxQR[iRight].y;
			ptBottom.x = ptCtrBoxQR[iBottom].x;
			ptTop.x = ptCtrBoxQR[iTop].x;

			double dTheta1 = asin(abs(ptLeft.y - ptCenter.y) / dRadius);
			ptLeft.x = ptCenter.x - dRadius * cos(dTheta1);
			ptLeft1.x -= iXDiff;

			double dTheta2 = asin(abs(ptRight.y - ptCenter.y) / dRadius);
			ptRight.x = ptCenter.x + dRadius * cos(dTheta2);
			ptRight1.x += iXDiff;

			double dTheta3 = acos(abs(ptBottom.x - ptCenter.x) / dRadius);
			ptBottom.y = ptCenter.y + dRadius * sin(dTheta3);
			ptBottom1.y += iYDiff;

			double dTheta4 = acos(abs(ptTop.x - ptCenter.x) / dRadius);
			ptTop.y = ptCenter.y - dRadius * sin(dTheta3);
			ptTop1.y -= iYDiff;

			iTopDistance = abs(ptCtrBoxQR[iTop].y - ptTop.y);
			iLeftDistance = abs(ptCtrBoxQR[iLeft].x - ptLeft.x);
			iRightDistance = abs(ptCtrBoxQR[iRight].x - ptRight.x);
			iBottomDistance = abs(ptCtrBoxQR[iBottom].y - ptBottom.y);

			line(originalImage, ptCtrBoxQR[iLeft], ptLeft, cvScalar(220, 0, 0), 2);
			line(originalImage, ptCtrBoxQR[iRight], ptRight, cvScalar(220, 0, 0), 2);
			line(originalImage, ptCtrBoxQR[iBottom], ptBottom, cvScalar(220, 0, 0), 2);
			line(originalImage, ptCtrBoxQR[iTop], ptTop, cvScalar(220, 0, 0), 2);
		}

		imwrite(strOutputImagePath, originalImage);
		////m_ofs << rstrImageFile << "," << (dRadius*PIXEL_TO_MM) << "," << (iLeftDistance*PIXEL_TO_MM) << "," << (iTopDistance*PIXEL_TO_MM) << "," << (iRightDistance*PIXEL_TO_MM) << "," << (iBottomDistance*PIXEL_TO_MM) << std::endl;
		m_ofs << rstrImageFile << "," << dRadius << "," << iLeftDistance << "," << iTopDistance << "," << iRightDistance << "," << iBottomDistance << std::endl;


		bRetVal = true;

		break;

	}

	if (storage) {
		cvReleaseMemStorage(&storage);
		storage = nullptr;
	}
	if (imgGrayScale) {
		cvReleaseImage(&imgGrayScale);
		imgGrayScale = nullptr;
	}


	return bRetVal;
}



bool QRCodeDetector::Process(const std::string& rstrImageFile)
{
	bool			bRetVal			= false;
	IplImage*		img				= nullptr;
	IplImage*		imgGrayScale	= nullptr;
	CvMemStorage*	storage			= nullptr;
	std::string		strImagePath{ m_strRootDirectory + rstrImageFile};
	std::string		strOutputImagePath{ m_strOutputDirectory + rstrImageFile };

	while(true) {

		cv::Mat image = cv::imread(strImagePath, 0); // IMG_3421.jpeg", 0);sample.jpg
		if (image.empty()) {
			break;
		}

		img = cvLoadImage(strImagePath.c_str());


		GaussianBlur(image, image, cv::Size(9, 9), 2, 2);

		//	//medianBlur(img, img, 5);
		//
		double dRadius = -1;
		std::vector<cv::Vec3f> circles;
		cv::Point ptCenter;
		ptCenter.x = -1;
		ptCenter.y = -1;

		////HoughCircles(image, circles, cv::HOUGH_GRADIENT, 1, image.rows/ 4, 200, 100); // change the last two parameters
		////HoughCircles(image, circles, cv::HOUGH_GRADIENT, 1.2, 200); // change the last two parameters
		////HoughCircles(image, circles, cv::HOUGH_GRADIENT, 1, 200); // change the last two parameters
		////HoughCircles(image, circles, cv::HOUGH_GRADIENT, 1, 250); // change the last two parameters
		////HoughCircles(image, circles, cv::HOUGH_GRADIENT, 1, 300); // change the last two parameters
		////HoughCircles(image, circles, cv::HOUGH_GRADIENT, 1, 100, 200, 100); // change the last two parameters

		HoughCircles(image, circles, cv::HOUGH_GRADIENT, 1, image.rows / 4, 200, 100, 0, 800); // change the last two parameters



		if (circles.size() <= 0) {
			cvSaveImage(strOutputImagePath.c_str(), img);
			m_ofs << rstrImageFile << "," << 0 << "," << 0 << "," << 0 << "," << 0 << std::endl;
			break;
		}

		for (size_t i = 0; i < circles.size(); i++) {
			cv::Vec3i c = circles[i];
			cvCircle(img, cv::Point(c[0], c[1]), c[2],cv::Scalar(0, 0, 255), 3, cv::LINE_AA);
			cvCircle(img, cv::Point(c[0], c[1]), 2, cv::Scalar(255, 0, 0), 3, cv::LINE_AA);

			dRadius = c[2];
			ptCenter.x = c[0];
			ptCenter.y = c[1];

			/*Point center(cvRound(circles[i][0]), cvRound(circles[i][1]));
			int radius = cvRound(circles[i][2]);
			circle(image, center, 3, Scalar(0, 255, 0), -1, 8, 0);
			circle(image, center, radius, Scalar(0, 0, 255), 3, 8, 0);
			cv::Vec2f v1= circles[i][0];
			cv::Vec2f v2 = circles[i][1];
			cv::Vec2f v3 = circles[i][2];

			reader.possibleCenters.push_back(v1);
			reader.possibleCenters.push_back(v2);
			reader.possibleCenters.push_back(v2);*/

		}

		//converting the original image into grayscale
		imgGrayScale = cvCreateImage(cvGetSize(img), 8, 1);
		cvCvtColor(img, imgGrayScale, CV_BGR2GRAY);

		//thresholding the grayscale image to get better results
		cvThreshold(imgGrayScale, imgGrayScale, 128, 255, CV_THRESH_BINARY);

		//converting the original image into grayscale
		//IplImage imgGrayScale(im); // cvCreateImage(cvGetSize(&crop), 8, 1);

		//thresholding the grayscale image to get better results
		//cvThreshold(&imgGrayScale, &imgGrayScale, 128, 255, CV_THRESH_BINARY);
		CvSeq*			contours = nullptr;  //hold the pointer to a contour in the memory block
		CvSeq*			result	= nullptr;   //hold sequence of points of a contour
		CvMemStorage*	storage = cvCreateMemStorage(0); //storage area for all contours

													   //finding all contours in the image
		cvFindContours(imgGrayScale, storage, &contours, sizeof(CvContour), CV_RETR_LIST, CV_CHAIN_APPROX_SIMPLE, cvPoint(0, 0));
		cv::Point ptCtrBoxQR[10];
		double dArmOR[10];
		int iQR = 0;
		int iYDiff = 0;
		int iXDiff = 0;
		int iTopDistance = 0;
		int iLeftDistance = 0;
		int iRightDistance = 0;
		int iBottomDistance = 0;
		double dArm = 0;
		double dArmMax = 0;
		double dArmMin = 9999;

		//iterating through each contour
		while (contours)
		{			//obtain a sequence of points of contour, pointed by the variable 'contour'
			result = cvApproxPoly(contours, sizeof(CvContour), storage, CV_POLY_APPROX_DP, cvContourPerimeter(contours)*0.02, 0);

			//if there are 3  vertices  in the contour(It should be a triangle)
			//if (result->total == 3)
			//{
			//	//iterating through each point
			//	CvPoint *pt[3];
			//	for (int i = 0; i < 3; i++) {
			//		pt[i] = (CvPoint*)cvGetSeqElem(result, i);
			//	}

			//	//drawing lines around the triangle
			//	cvLine(img, *pt[0], *pt[1], cvScalar(255, 0, 0), 4);
			//	cvLine(img, *pt[1], *pt[2], cvScalar(255, 0, 0), 4);
			//	cvLine(img, *pt[2], *pt[0], cvScalar(255, 0, 0), 4);
			//}

			//if there are 4 vertices in the contour(It should be a quadrilateral)
			CvPoint *pt[4];
			cv::Point ptCtrBox;
			ptCtrBox.x = -1;
			ptCtrBox.y = -1;


			if (result->total == 4)
			{
				//iterating through each point
				for (int i = 0; i < 4; i++) {
					pt[i] = (CvPoint*)cvGetSeqElem(result, i);
				}

				double dPtDistance = sqrt(
										(((*pt[0]).x - ptCenter.x)* ((*pt[0]).x - ptCenter.x)) +
										(((*pt[0]).y - ptCenter.y)* ((*pt[0]).y - ptCenter.y))
										);
				if (dPtDistance < dRadius) {
					if ((((abs((*pt[0]).x - (*pt[1]).x) > 15) ||
						(abs((*pt[0]).y - (*pt[1]).y) > 15)) &&
						((abs((*pt[1]).x - (*pt[2]).x) > 15) ||
						(abs((*pt[1]).y - (*pt[2]).y) > 15)) &&
							((abs((*pt[2]).x - (*pt[3]).x) > 15) ||
						(abs((*pt[2]).y - (*pt[3]).y) > 15)) &&
								((abs((*pt[3]).x - (*pt[0]).x) > 15) ||
						(abs((*pt[3]).y - (*pt[0]).y) > 15))) &&

									(((abs((*pt[0]).x - (*pt[1]).x) < 60) ||
						(abs((*pt[0]).y - (*pt[1]).y) < 60)) &&
										((abs((*pt[1]).x - (*pt[2]).x) < 60) ||
										(abs((*pt[1]).y - (*pt[2]).y) < 60)) &&
											((abs((*pt[2]).x - (*pt[3]).x) < 60) ||
										(abs((*pt[2]).y - (*pt[3]).y) < 60)) &&
												((abs((*pt[3]).x - (*pt[0]).x) < 60) ||
										(abs((*pt[3]).y - (*pt[0]).y) < 60))))
					{

						double dArm1 = sqrt(
							(((*pt[0]).x - (*pt[1]).x)* ((*pt[0]).x - (*pt[1]).x)) +
							(((*pt[0]).y - (*pt[1]).y)* ((*pt[0]).y - (*pt[1]).y))
							);
						double dArm2 = sqrt(
							(((*pt[1]).x - (*pt[2]).x)* ((*pt[1]).x - (*pt[2]).x)) +
							(((*pt[1]).y - (*pt[2]).y)* ((*pt[1]).y - (*pt[2]).y))
						);
						double dArm3 = sqrt(
							(((*pt[2]).x - (*pt[3]).x)* ((*pt[2]).x - (*pt[3]).x)) +
							(((*pt[2]).y - (*pt[3]).y)* ((*pt[2]).y - (*pt[3]).y))
						);
						double dArm4 = sqrt(
							(((*pt[3]).x - (*pt[0]).x)* ((*pt[3]).x - (*pt[0]).x)) +
							(((*pt[3]).y - (*pt[0]).y)* ((*pt[3]).y - (*pt[0]).y))
						);
						dArm = (dArm1 + dArm2 + dArm3 + dArm4)/4;
						if (dArmMax < dArm)
							dArmMax = dArm;
						if (dArmMin > dArm)
							dArmMin = dArm;
						cvLine(img, *pt[0], *pt[1], cvScalar(0, 255, 0), 4);
						cvLine(img, *pt[1], *pt[2], cvScalar(0, 255, 0), 4);
						cvLine(img, *pt[2], *pt[3], cvScalar(0, 255, 0), 4);
						cvLine(img, *pt[3], *pt[0], cvScalar(0, 255, 0), 4);

						ptCtrBox.x = double((*pt[0]).x + (*pt[1]).x + (*pt[2]).x + (*pt[3]).x) / 4;
						ptCtrBox.y = double((*pt[0]).y + (*pt[1]).y + (*pt[2]).y + (*pt[3]).y) / 4;

						iYDiff = abs(ptCtrBox.y - (*pt[0]).y);
						iXDiff = abs(ptCtrBox.x - (*pt[0]).x);

						bool bExist = false;
						if (iQR > 0) {
							for (int iQ = 0; iQ < iQR; iQ++) {
								if ((abs(ptCtrBoxQR[iQ].x - ptCtrBox.x) < 5) && (abs(ptCtrBoxQR[iQ].y - ptCtrBox.y) < 5)) {
									bExist = true;
								}
							}
						}

						if (!bExist) {
							ptCtrBoxQR[iQR].x = ptCtrBox.x;
							ptCtrBoxQR[iQR].y = ptCtrBox.y;
							dArmOR[iQR] = dArm;
							iQR = iQR + 1;
						}
					}
				}
			}

			//if there are 7  vertices  in the contour(It should be a heptagon)
			//else if (result->total == 7)
			//{
			//	//iterating through each point
			//	CvPoint *pt[7];
			//	for (int i = 0; i < 7; i++) {
			//		pt[i] = (CvPoint*)cvGetSeqElem(result, i);
			//	}
			//	cvLine(img, *pt[0], *pt[1], cvScalar(0, 0, 255), 4);
			//	cvLine(img, *pt[1], *pt[2], cvScalar(0, 0, 255), 4);
			//	cvLine(img, *pt[2], *pt[3], cvScalar(0, 0, 255), 4);
			//	cvLine(img, *pt[3], *pt[4], cvScalar(0, 0, 255), 4);
			//	cvLine(img, *pt[4], *pt[5], cvScalar(0, 0, 255), 4);
			//	cvLine(img, *pt[5], *pt[6], cvScalar(0, 0, 255), 4);
			//	cvLine(img, *pt[6], *pt[0], cvScalar(0, 0, 255), 4);
			//	//drawing lines around the heptagon
			//	//cv::line(im, *pt[0], *pt[1], cvScalar(0, 0, 255), 4);
			//	//cv::line(im, *pt[1], *pt[2], cvScalar(0, 0, 255), 4);
			//	//cv::line(im, *pt[2], *pt[3], cvScalar(0, 0, 255), 4);
			//	//cv::line(im, *pt[3], *pt[4], cvScalar(0, 0, 255), 4);
			//	//cv::line(im, *pt[4], *pt[5], cvScalar(0, 0, 255), 4);
			//	//cv::line(im, *pt[5], *pt[6], cvScalar(0, 0, 255), 4);
			//	//cv::line(im, *pt[6], *pt[0], cvScalar(0, 0, 255), 4);
			//}
			//obtain the next contour

			if ((dRadius != -1) && (ptCtrBox.x != -1) && (ptCtrBox.y != -1)) {
				double dDistanceBoxFromCtr = sqrt(((ptCtrBox.x - ptCenter.x) * (ptCtrBox.x - ptCenter.x)) +
					((ptCtrBox.y - ptCenter.y) * (ptCtrBox.y - ptCenter.y)));

				double dDistanceFromOutside = dRadius - dDistanceBoxFromCtr;

				cvLine(img, ptCtrBox, ptCenter, cvScalar(0, 255, 0), 2);

				//	cv::waitKey(0);
				int i = 0;
				i++;
			}

			contours = contours->h_next;
		}

		if (iQR >= 3) {
			if (iQR > 3) { // elemenate south west square
				while (true) {
					int iS1 = -1;
					int iS2 = -1;
					double dMax1 = 0;

					//if((ptCtrBoxQR[0].y > ptCtrBoxQR[1].y)
					for (int iE = 0; iE < iQR; iE++) {
						if (dMax1 < ptCtrBoxQR[iE].y) {
							dMax1 = ptCtrBoxQR[iE].y;
							iS1 = iE;
						}
					}

					dMax1 = 0;

					//if((ptCtrBoxQR[0].y > ptCtrBoxQR[1].y)
					for (int iE = 0; iE < iQR; iE++) {
						if (iE != iS1) {
							if (dMax1 < ptCtrBoxQR[iE].y) {
								dMax1 = ptCtrBoxQR[iE].y;
								iS2 = iE;
							}
						}
					}

					int iSW = -1;
					if ((iS1 >= 0) && (iS2 >= 0)) {
						if (ptCtrBoxQR[iS1].x > ptCtrBoxQR[iS2].x) {
							iSW = iS1;
						}
						else {
							iSW = iS2;
						}
					}

					if (iSW >= 0) {
						ptCtrBoxQR[iSW].x = -1;
						ptCtrBoxQR[iSW].y = -1;
						for (int iE = iSW; iE < iQR; iE++) {
							ptCtrBoxQR[iE].x = ptCtrBoxQR[iE + 1].x;
							ptCtrBoxQR[iE].y = ptCtrBoxQR[iE + 1].y;
							if (iE + 1 < 10) {
								ptCtrBoxQR[iE + 1].x = -1;
								ptCtrBoxQR[iE + 1].y = -1;
							}
						}
						iQR--;
					}
					if (iQR == 3) {
						break;
					}
				}
			}
			
			/*for (int iE = 0; iE < iQR - 1; iE++) {
				if ((dArmOR[iE] == 0) && (ptCtrBoxQR[iE].x == 0) && (ptCtrBoxQR[iE].y == 0)) {
					break;
				}
				if (dArmOR[iE] == dArmMin) {
					dArmOR[iE] = dArmOR[iE + 1];
					ptCtrBoxQR[iE].x = ptCtrBoxQR[iE + 1].x;
					ptCtrBoxQR[iE].y = ptCtrBoxQR[iE + 1].y;
				}
			}
			while (true) {
				int iElem = -1;
				for (int iE = 0; iE < iQR - 1; iE++) {
					if ((dArmOR[iE] == 0) && (ptCtrBoxQR[iE].x == 0) && (ptCtrBoxQR[iE].y == 0)) {
						break;
					}
					if (abs(dArmMax - dArmOR[iQR]) < 5) {
					}
					else {
						iElem = iE;
						dArmOR[iE] = dArmOR[iE + 1];
						ptCtrBoxQR[iE].x = ptCtrBoxQR[iE + 1].x;
						ptCtrBoxQR[iE].y = ptCtrBoxQR[iE + 1].y;
					}
				}
				if (iElem == -1) {
					break;
				}
			}*/

			int iBottom = -1;
			int iLeft = -1;
			int iRight = -1;
			int iTop = -1;

			cv::Point ptLeft1;
			cv::Point ptTop1;
			cv::Point ptRight1;
			cv::Point ptBottom1;

			ptLeft1.x = -1;
			ptLeft1.y = -1;

			ptTop1.x = -1;
			ptTop1.y = -1;

			ptRight1.x = -1;
			ptRight1.y = -1;

			ptBottom1.x = -1;
			ptBottom1.y = -1;

			if ((abs(ptCtrBoxQR[2].y - ptCtrBoxQR[0].y) > dArm * 3) && (abs(ptCtrBoxQR[2].y - ptCtrBoxQR[1].y) > dArm * 3) &&
				(abs(ptCtrBoxQR[1].y - ptCtrBoxQR[0].y) < dArm * 2)) {
				iBottom = 2;

				if (ptCtrBoxQR[1].x > ptCtrBoxQR[0].x) {
					iLeft = 0;
					iRight = 1;
				}
				else {
					iLeft = 1;
					iRight = 0;
				}
				if (ptCtrBoxQR[1].y > ptCtrBoxQR[0].y) {
					iTop = 0;
				}
				else {
					iTop = 1;
				}

			}
			else if ((abs(ptCtrBoxQR[1].y - ptCtrBoxQR[0].y) > dArm * 3) && (abs(ptCtrBoxQR[1].y - ptCtrBoxQR[1].y) > dArm * 3) &&
				(abs(ptCtrBoxQR[2].y - ptCtrBoxQR[0].y) < dArm * 2)) {
				iBottom = 1;
				if (ptCtrBoxQR[2].x > ptCtrBoxQR[0].x) {
					iLeft = 0;
					iRight = 2;
				}
				else {
					iLeft = 2;
					iRight = 0;
				}
				if (ptCtrBoxQR[2].y > ptCtrBoxQR[0].y) {
					iTop = 0;
				}
				else {
					iTop = 1;
				}

			}
			else if ((abs(ptCtrBoxQR[0].y - ptCtrBoxQR[1].y) > dArm * 3) && (abs(ptCtrBoxQR[0].y - ptCtrBoxQR[2].y) > dArm * 3) &&
				(abs(ptCtrBoxQR[2].y - ptCtrBoxQR[1].y) < dArm * 2)) {
				iBottom = 0;
				if (ptCtrBoxQR[2].x > ptCtrBoxQR[1].x) {
					iLeft = 1;
					iRight = 2;
				}
				else {
					iLeft = 2;
					iRight = 1;
				}
				if (ptCtrBoxQR[1].y > ptCtrBoxQR[1].y) {
					iTop = 1;
				}
				else {
					iTop = 2;
				}

			}
			if (iLeft > -1 && iRight > -1 && iBottom > -1 && iTop>-1) {
				if ((abs(ptCtrBoxQR[iLeft].y - ptCtrBoxQR[iRight].y)) < 30) {
					ptTop1.x = (ptCtrBoxQR[iLeft].x + ptCtrBoxQR[iRight].x) / 2;
					ptTop1.y = (ptCtrBoxQR[iLeft].y + ptCtrBoxQR[iRight].y) / 2;
				}
				else if (ptCtrBoxQR[iLeft].y < ptCtrBoxQR[iRight].y) {
					ptTop1.x = ptCtrBoxQR[iLeft].x;
					ptTop1.y = ptCtrBoxQR[iLeft].y;
				}
				else {
					ptTop1.x = ptCtrBoxQR[iRight].x;
					ptTop1.y = ptCtrBoxQR[iRight].y;
				}

				if ((abs(ptCtrBoxQR[iLeft].x - ptCtrBoxQR[iBottom].x)) < 30) {
					ptLeft1.x = (ptCtrBoxQR[iLeft].x + ptCtrBoxQR[iBottom].x) / 2;
					ptLeft1.y = (ptCtrBoxQR[iLeft].y + ptCtrBoxQR[iBottom].y) / 2;
				}
				else if (ptCtrBoxQR[iLeft].x < ptCtrBoxQR[iBottom].x) {
					ptLeft1.x = ptCtrBoxQR[iLeft].x;
					ptLeft1.y = ptCtrBoxQR[iLeft].y;
				}

				if ((abs(ptCtrBoxQR[iRight].x - ptCtrBoxQR[iBottom].x)) < 30) {
					ptRight1.x = (ptCtrBoxQR[iRight].x + ptCtrBoxQR[iBottom].x) / 2;
					ptRight1.y = (ptCtrBoxQR[iRight].y + ptCtrBoxQR[iBottom].y) / 2;
				}
				else if (ptCtrBoxQR[iRight].x > ptCtrBoxQR[iBottom].x) {
					ptRight1.x = ptCtrBoxQR[iRight].x;
					ptRight1.y = ptCtrBoxQR[iRight].y;
				}
				/*else {
				ptRight.x = ptCtrBoxQR[2].x;
				ptRight.y = ptCtrBoxQR[2].y;
				}*/

				if ((abs(ptCtrBoxQR[iBottom].y - ptCtrBoxQR[iRight].y)) < 30) {
					ptBottom1.x = (ptCtrBoxQR[iBottom].x + ptCtrBoxQR[iRight].x) / 2;
					ptBottom1.y = (ptCtrBoxQR[iBottom].y + ptCtrBoxQR[iRight].y) / 2;
				}
				else if ((abs(ptCtrBoxQR[iBottom].y - ptCtrBoxQR[iLeft].y)) < 30) {
					ptBottom1.x = (ptCtrBoxQR[iBottom].x + ptCtrBoxQR[iLeft].x) / 2;
					ptBottom1.y = (ptCtrBoxQR[iBottom].y + ptCtrBoxQR[iLeft].y) / 2;
				}
				else {
					ptBottom1.x = ptCtrBoxQR[iBottom].x;
					ptBottom1.y = ptCtrBoxQR[iBottom].y;
				}
			}

			cv::Point ptLeft;
			cv::Point ptRight;
			cv::Point ptBottom;
			cv::Point ptTop;

			//ptLeft.y = ptLeft1.y;
			//ptRight.y = ptRight1.y;
			//ptBottom.x = ptBottom1.x;
			//ptTop.x = ptTop1.x;

			ptLeft.y = ptCtrBoxQR[iLeft].y;
			ptRight.y = ptCtrBoxQR[iRight].y;
			ptBottom.x = ptCtrBoxQR[iBottom].x;
			ptTop.x = ptCtrBoxQR[iTop].x;

			double dTheta1 = asin(abs(ptLeft.y - ptCenter.y) / dRadius);
			ptLeft.x = ptCenter.x - dRadius * cos(dTheta1);
			ptLeft1.x -= iXDiff;

			double dTheta2 = asin(abs(ptRight.y - ptCenter.y) / dRadius);
			ptRight.x = ptCenter.x + dRadius * cos(dTheta2);
			ptRight1.x += iXDiff;

			double dTheta3 = acos(abs(ptBottom.x - ptCenter.x) / dRadius);
			ptBottom.y = ptCenter.y + dRadius * sin(dTheta3);
			ptBottom1.y += iYDiff;

			double dTheta4 = acos(abs(ptTop.x - ptCenter.x) / dRadius);
			ptTop.y = ptCenter.y - dRadius * sin(dTheta3);
			ptTop1.y -= iYDiff;

			iTopDistance = abs(ptCtrBoxQR[iTop].y - ptTop.y);
			iLeftDistance = abs(ptCtrBoxQR[iLeft].x - ptLeft.x);
			iRightDistance = abs(ptCtrBoxQR[iRight].x - ptRight.x);
			iBottomDistance = abs(ptCtrBoxQR[iBottom].y - ptBottom.y);

			cvLine(img, ptCtrBoxQR[iLeft], ptLeft, cvScalar(220, 0, 0), 2);
			cvLine(img, ptCtrBoxQR[iRight], ptRight, cvScalar(220, 0, 0), 2);
			cvLine(img, ptCtrBoxQR[iBottom], ptBottom, cvScalar(220, 0, 0), 2);
			cvLine(img, ptCtrBoxQR[iTop], ptTop, cvScalar(220, 0, 0), 2);

			//iTopDistance = abs(ptTop1.y - ptTop.y);
			//iLeftDistance = abs(ptLeft1.x - ptLeft.x);
			//iRightDistance = abs(ptRight1.x - ptRight.x);
			//iBottomDistance = abs(ptBottom1.y - ptBottom.y);

			//cvLine(img, ptLeft1, ptLeft, cvScalar(220, 0, 0), 2);
			//cvLine(img, ptRight1, ptRight, cvScalar(220, 0, 0), 2);
			//cvLine(img, ptBottom1, ptBottom, cvScalar(220, 0, 0), 2);
			//cvLine(img, ptTop1, ptTop, cvScalar(220, 0, 0), 2);

			//if (((ptCtrBoxQR[2].y - ptCtrBoxQR[0].y) > 50) && ((ptCtrBoxQR[2].y - ptCtrBoxQR[1].y) > 50) &&
			//	((ptCtrBoxQR[1].y - ptCtrBoxQR[0].y) < 30)) {
			//	bTop2 = true;
			//}

			////QRDistance QrDist;
			////QrDist.iLeft_	= iLeftDistance;
			////QrDist.iRight_	= iRightDistance;
			////QrDist.iTop_	= iTopDistance;
			////QrDist.iBottom_ = iBottomDistance;

			////if (pWndParent_) {
			////	if (pWndParent_->GetSafeHwnd()) {
			////		pWndParent_->SendMessage(1001, 0, LPARAM(&QrDist));
			////	}
			////}
		}
		cvNamedWindow(strImagePath.c_str(), CV_WINDOW_NORMAL);
		cvShowImage(strImagePath.c_str(), img);

		cvSaveImage(strOutputImagePath.c_str(), img);
		m_ofs << rstrImageFile << "," << (iLeftDistance*PIXEL_TO_MM) << "," << (iTopDistance*PIXEL_TO_MM) << "," << (iRightDistance*PIXEL_TO_MM) << "," << (iBottomDistance*PIXEL_TO_MM) << std::endl;

		bRetVal = true;

		break;

	}

	if (storage) {
		cvReleaseMemStorage(&storage);
		storage = nullptr;
	}
	if (imgGrayScale) {
		cvReleaseImage(&imgGrayScale);
		imgGrayScale = nullptr;
	}

	if (img) {
		cvReleaseImage(&img);
		img = nullptr;
	}

	cvDestroyWindow(strImagePath.c_str());

	return bRetVal;
}


////bool QRCodeDetector::Process()
////{
////
////	cv::Mat image = cv::imread(m_strImagePath, 0);
////	IplImage* img = cvLoadImage(m_strImagePath.c_str());
////
////	//show the original image
////	//cvNamedWindow("Raw");
////	//cvShowImage("Raw", img);
////
////	//cvtColor(image, image, CV_BGR2GRAY);
////
////	GaussianBlur(image, image, cv::Size(11, 11), 2, 2);
////
////	//	//medianBlur(img, img, 5);
////	//
////	double dRadius = -1;
////	std::vector<cv::Vec3f> circles;
////	cv::Point ptCenter;
////	ptCenter.x = -1;
////	ptCenter.y = -1;
////
////	HoughCircles(image, circles, cv::HOUGH_GRADIENT, 1, image.rows / 4, 200, 100); // change the last two parameters
////
////	if (circles.size() <= 0) {
////		return false;
////	}
////
////	for (size_t i = 0; i < circles.size(); i++) {
////		cv::Vec3i c = circles[i];
////		cvCircle(img, cv::Point(c[0], c[1]), c[2], cv::Scalar(0, 0, 255), 3, cv::LINE_AA);
////		cvCircle(img, cv::Point(c[0], c[1]), 2, cv::Scalar(255, 0, 0), 3, cv::LINE_AA);
////
////		dRadius = c[2];
////		ptCenter.x = c[0];
////		ptCenter.y = c[1];
////
////		/*Point center(cvRound(circles[i][0]), cvRound(circles[i][1]));
////		int radius = cvRound(circles[i][2]);
////		circle(image, center, 3, Scalar(0, 255, 0), -1, 8, 0);
////		circle(image, center, radius, Scalar(0, 0, 255), 3, 8, 0);
////		cv::Vec2f v1= circles[i][0];
////		cv::Vec2f v2 = circles[i][1];
////		cv::Vec2f v3 = circles[i][2];
////
////		reader.possibleCenters.push_back(v1);
////		reader.possibleCenters.push_back(v2);
////		reader.possibleCenters.push_back(v2);*/
////
////	}
////
////	//converting the original image into grayscale
////	IplImage* imgGrayScale = cvCreateImage(cvGetSize(img), 8, 1);
////	cvCvtColor(img, imgGrayScale, CV_BGR2GRAY);
////
////	//thresholding the grayscale image to get better results
////	cvThreshold(imgGrayScale, imgGrayScale, 128, 255, CV_THRESH_BINARY);
////
////	//converting the original image into grayscale
////	//IplImage imgGrayScale(im); // cvCreateImage(cvGetSize(&crop), 8, 1);
////
////	//thresholding the grayscale image to get better results
////	//cvThreshold(&imgGrayScale, &imgGrayScale, 128, 255, CV_THRESH_BINARY);
////	CvSeq* contours = nullptr;  //hold the pointer to a contour in the memory block
////	CvSeq* result	= nullptr;   //hold sequence of points of a contour
////	CvMemStorage *storage = cvCreateMemStorage(0); //storage area for all contours
////
////												   //finding all contours in the image
////	cvFindContours(imgGrayScale, storage, &contours, sizeof(CvContour), CV_RETR_LIST, CV_CHAIN_APPROX_SIMPLE, cvPoint(0, 0));
////	cv::Point ptCtrBoxQR[30];
////	int iQR = 0;
////
////	//iterating through each contour
////	while (contours)
////	{
////		//obtain a sequence of points of contour, pointed by the variable 'contour'
////		result = cvApproxPoly(contours, sizeof(CvContour), storage, CV_POLY_APPROX_DP, cvContourPerimeter(contours)*0.02, 0);
////
////		//if there are 3  vertices  in the contour(It should be a triangle)
////		//if (result->total == 3)
////		//{
////		//	//iterating through each point
////		//	CvPoint *pt[3];
////		//	for (int i = 0; i < 3; i++) {
////		//		pt[i] = (CvPoint*)cvGetSeqElem(result, i);
////		//	}
////
////		//	//drawing lines around the triangle
////		//	cvLine(img, *pt[0], *pt[1], cvScalar(255, 0, 0), 4);
////		//	cvLine(img, *pt[1], *pt[2], cvScalar(255, 0, 0), 4);
////		//	cvLine(img, *pt[2], *pt[0], cvScalar(255, 0, 0), 4);
////		//}
////
////		//if there are 4 vertices in the contour(It should be a quadrilateral)
////		CvPoint *pt[4];
////		cv::Point ptCtrBox;
////		ptCtrBox.x = -1;
////		ptCtrBox.y = -1;
////
////		if (result->total == 4)
////		{
////			//iterating through each point
////			for (int i = 0; i < 4; i++) {
////				pt[i] = (CvPoint*)cvGetSeqElem(result, i);
////			}
////			if ((((abs((*pt[0]).x - (*pt[1]).x) > 15) ||
////				(abs((*pt[0]).y - (*pt[1]).y) > 15)) &&
////				((abs((*pt[1]).x - (*pt[2]).x) > 15) ||
////				(abs((*pt[1]).y - (*pt[2]).y) > 15)) &&
////					((abs((*pt[2]).x - (*pt[3]).x) > 15) ||
////				(abs((*pt[2]).y - (*pt[3]).y) > 15)) &&
////						((abs((*pt[3]).x - (*pt[0]).x) > 5) ||
////				(abs((*pt[3]).y - (*pt[0]).y) > 5))) &&
////
////							(((abs((*pt[0]).x - (*pt[1]).x) < 30) ||
////				(abs((*pt[0]).y - (*pt[1]).y) < 30)) &&
////								((abs((*pt[1]).x - (*pt[2]).x) < 30) ||
////								(abs((*pt[1]).y - (*pt[2]).y) < 30)) &&
////									((abs((*pt[2]).x - (*pt[3]).x) < 30) ||
////								(abs((*pt[2]).y - (*pt[3]).y) < 30)) &&
////										((abs((*pt[3]).x - (*pt[0]).x) < 30) ||
////								(abs((*pt[3]).y - (*pt[0]).y) < 30))))
////			{
////
////				cvLine(img, *pt[0], *pt[1], cvScalar(0, 255, 0), 4);
////				cvLine(img, *pt[1], *pt[2], cvScalar(0, 255, 0), 4);
////				cvLine(img, *pt[2], *pt[3], cvScalar(0, 255, 0), 4);
////				cvLine(img, *pt[3], *pt[0], cvScalar(0, 255, 0), 4);
////
////				ptCtrBox.x = double((*pt[0]).x + (*pt[1]).x + (*pt[2]).x + (*pt[3]).x) / 4;
////				ptCtrBox.y = double((*pt[0]).y + (*pt[1]).y + (*pt[2]).y + (*pt[3]).y) / 4;
////
////				ptCtrBoxQR[iQR].x = ptCtrBox.x;
////				ptCtrBoxQR[iQR++].y = ptCtrBox.y;
////			}
////		}
////
////		//if there are 7  vertices  in the contour(It should be a heptagon)
////		//else if (result->total == 7)
////		//{
////		//	//iterating through each point
////		//	CvPoint *pt[7];
////		//	for (int i = 0; i < 7; i++) {
////		//		pt[i] = (CvPoint*)cvGetSeqElem(result, i);
////		//	}
////		//	cvLine(img, *pt[0], *pt[1], cvScalar(0, 0, 255), 4);
////		//	cvLine(img, *pt[1], *pt[2], cvScalar(0, 0, 255), 4);
////		//	cvLine(img, *pt[2], *pt[3], cvScalar(0, 0, 255), 4);
////		//	cvLine(img, *pt[3], *pt[4], cvScalar(0, 0, 255), 4);
////		//	cvLine(img, *pt[4], *pt[5], cvScalar(0, 0, 255), 4);
////		//	cvLine(img, *pt[5], *pt[6], cvScalar(0, 0, 255), 4);
////		//	cvLine(img, *pt[6], *pt[0], cvScalar(0, 0, 255), 4);
////		//	//drawing lines around the heptagon
////		//	//cv::line(im, *pt[0], *pt[1], cvScalar(0, 0, 255), 4);
////		//	//cv::line(im, *pt[1], *pt[2], cvScalar(0, 0, 255), 4);
////		//	//cv::line(im, *pt[2], *pt[3], cvScalar(0, 0, 255), 4);
////		//	//cv::line(im, *pt[3], *pt[4], cvScalar(0, 0, 255), 4);
////		//	//cv::line(im, *pt[4], *pt[5], cvScalar(0, 0, 255), 4);
////		//	//cv::line(im, *pt[5], *pt[6], cvScalar(0, 0, 255), 4);
////		//	//cv::line(im, *pt[6], *pt[0], cvScalar(0, 0, 255), 4);
////		//}
////		//obtain the next contour
////
////		if ((dRadius != -1) && (ptCtrBox.x != -1) && (ptCtrBox.y != -1)) {
////			double dDistanceBoxFromCtr = sqrt(((ptCtrBox.x - ptCenter.x) * (ptCtrBox.x - ptCenter.x)) +
////				((ptCtrBox.y - ptCenter.y) * (ptCtrBox.y - ptCenter.y)));
////
////			double dDistanceFromOutside = dRadius - dDistanceBoxFromCtr;
////
////			cvLine(img, ptCtrBox, ptCenter, cvScalar(0, 255, 0), 2);
////
////			//	cv::waitKey(0);
////			int i = 0;
////			i++;
////		}
////
////		contours = contours->h_next;
////	}
////
////	if (iQR >= 3) {
////		int iBottom = -1;
////		int iLeft = -1;
////		int iRight = -1;
////
////		cv::Point ptLeft1;
////		cv::Point ptTop1;
////		cv::Point ptRight1;
////		cv::Point ptBottom1;
////
////		ptLeft1.x = -1;
////		ptLeft1.y = -1;
////
////		ptTop1.x = -1;
////		ptTop1.y = -1;
////
////		ptRight1.x = -1;
////		ptRight1.y = -1;
////
////		ptBottom1.x = -1;
////		ptBottom1.y = -1;
////
////		if (((ptCtrBoxQR[2].y - ptCtrBoxQR[0].y) > 50) && ((ptCtrBoxQR[2].y - ptCtrBoxQR[1].y) > 50) &&
////			(abs(ptCtrBoxQR[1].y - ptCtrBoxQR[0].y) < 30)) {
////			iBottom = 2;
////
////			if (ptCtrBoxQR[1].x > ptCtrBoxQR[0].x) {
////				iLeft = 0;
////				iRight = 1;
////			}
////			else {
////				iLeft = 1;
////				iRight = 0;
////			}
////
////		}
////		else if (((ptCtrBoxQR[1].y - ptCtrBoxQR[0].y) > 50) && ((ptCtrBoxQR[1].y - ptCtrBoxQR[1].y) > 50) &&
////			(abs(ptCtrBoxQR[2].y - ptCtrBoxQR[0].y) < 30)) {
////			iBottom = 1;
////			if (ptCtrBoxQR[2].x > ptCtrBoxQR[0].x) {
////				iLeft = 0;
////				iRight = 2;
////			}
////			else {
////				iLeft = 2;
////				iRight = 0;
////			}
////		}
////		else if (((ptCtrBoxQR[0].y - ptCtrBoxQR[1].y) > 50) && ((ptCtrBoxQR[0].y - ptCtrBoxQR[2].y) > 50) &&
////			(abs(ptCtrBoxQR[2].y - ptCtrBoxQR[1].y) < 30)) {
////			iBottom = 0;
////			if (ptCtrBoxQR[2].x > ptCtrBoxQR[1].x) {
////				iLeft = 1;
////				iRight = 2;
////			}
////			else {
////				iLeft = 2;
////				iRight = 1;
////			}
////		}
////		if (iLeft > -1 && iRight > -1 && iBottom > -1) {
////			if ((abs(ptCtrBoxQR[iLeft].y - ptCtrBoxQR[iRight].y)) < 30) {
////				ptTop1.x = (ptCtrBoxQR[iLeft].x + ptCtrBoxQR[iRight].x) / 2;
////				ptTop1.y = (ptCtrBoxQR[iLeft].y + ptCtrBoxQR[iRight].y) / 2;
////			}
////			else if (ptCtrBoxQR[iLeft].y < ptCtrBoxQR[iRight].y) {
////				ptTop1.x = ptCtrBoxQR[iLeft].x;
////				ptTop1.y = ptCtrBoxQR[iLeft].y;
////			}
////			else {
////				ptTop1.x = ptCtrBoxQR[iRight].x;
////				ptTop1.y = ptCtrBoxQR[iRight].y;
////			}
////
////			if ((abs(ptCtrBoxQR[iLeft].x - ptCtrBoxQR[iBottom].x)) < 30) {
////				ptLeft1.x = (ptCtrBoxQR[iLeft].x + ptCtrBoxQR[iBottom].x) / 2;
////				ptLeft1.y = (ptCtrBoxQR[iLeft].y + ptCtrBoxQR[iBottom].y) / 2;
////			}
////			else if (ptCtrBoxQR[iLeft].x < ptCtrBoxQR[iBottom].x) {
////				ptLeft1.x = ptCtrBoxQR[iLeft].x;
////				ptLeft1.y = ptCtrBoxQR[iLeft].y;
////			}
////
////			if ((abs(ptCtrBoxQR[iRight].x - ptCtrBoxQR[iBottom].x)) < 30) {
////				ptRight1.x = (ptCtrBoxQR[iRight].x + ptCtrBoxQR[iBottom].x) / 2;
////				ptRight1.y = (ptCtrBoxQR[iRight].y + ptCtrBoxQR[iBottom].y) / 2;
////			}
////			else if (ptCtrBoxQR[iRight].x > ptCtrBoxQR[iBottom].x) {
////				ptRight1.x = ptCtrBoxQR[iRight].x;
////				ptRight1.y = ptCtrBoxQR[iRight].y;
////			}
////			/*else {
////			ptRight.x = ptCtrBoxQR[2].x;
////			ptRight.y = ptCtrBoxQR[2].y;
////			}*/
////
////			if ((abs(ptCtrBoxQR[iBottom].y - ptCtrBoxQR[iRight].y)) < 30) {
////				ptBottom1.x = (ptCtrBoxQR[iBottom].x + ptCtrBoxQR[iRight].x) / 2;
////				ptBottom1.y = (ptCtrBoxQR[iBottom].y + ptCtrBoxQR[iRight].y) / 2;
////			}
////			else if ((abs(ptCtrBoxQR[iBottom].y - ptCtrBoxQR[iLeft].y)) < 30) {
////				ptBottom1.x = (ptCtrBoxQR[iBottom].x + ptCtrBoxQR[iLeft].x) / 2;
////				ptBottom1.y = (ptCtrBoxQR[iBottom].y + ptCtrBoxQR[iLeft].y) / 2;
////			}
////			else {
////				ptBottom1.x = ptCtrBoxQR[iBottom].x;
////				ptBottom1.y = ptCtrBoxQR[iBottom].y;
////			}
////		}
////
////		cv::Point ptLeft;
////		cv::Point ptRight;
////		cv::Point ptBottom;
////		cv::Point ptTop;
////
////		ptLeft.y = ptLeft1.y;
////		ptRight.y = ptRight1.y;
////		ptBottom.x = ptBottom1.x;
////		ptTop.x = ptTop1.x;
////
////		double dTheta1 = asin(abs(ptLeft.y - ptCenter.y) / dRadius);
////		ptLeft.x = ptCenter.x - dRadius * cos(dTheta1);
////
////		double dTheta2 = asin(abs(ptRight.y - ptCenter.y) / dRadius);
////		ptRight.x = ptCenter.x + dRadius * cos(dTheta2);
////
////		double dTheta3 = acos(abs(ptBottom.x - ptCenter.x) / dRadius);
////		ptBottom.y = ptCenter.y + dRadius * sin(dTheta3);
////
////		double dTheta4 = acos(abs(ptTop.x - ptCenter.x) / dRadius);
////		ptTop.y = ptCenter.y - dRadius * sin(dTheta3);
////
////		//cvLine(img, ptCtrBoxQR[iLeft], ptLeft, cvScalar(0, 255, 0), 2);
////		//cvLine(img, ptCtrBoxQR[iRight], ptRight, cvScalar(0, 255, 0), 2);
////		//cvLine(img, ptCtrBoxQR[iBottom], ptBottom, cvScalar(0, 255, 0), 2);
////
////
////		cvLine(img, ptLeft1, ptLeft, cvScalar(220, 0, 0), 2);
////		cvLine(img, ptRight1, ptRight, cvScalar(220, 0, 0), 2);
////		cvLine(img, ptBottom1, ptBottom, cvScalar(220, 0, 0), 2);
////		cvLine(img, ptTop1, ptTop, cvScalar(220, 0, 0), 2);
////
////		//if (((ptCtrBoxQR[2].y - ptCtrBoxQR[0].y) > 50) && ((ptCtrBoxQR[2].y - ptCtrBoxQR[1].y) > 50) &&
////		//	((ptCtrBoxQR[1].y - ptCtrBoxQR[0].y) < 30)) {
////		//	bTop2 = true;
////		//}
////
////	}
////
////	//cv::imshow("QR code", im);
////	//cvShowImage("Tracked", img);
////	//cv::waitKey(10000);
////
////
////
////	cvNamedWindow("circles", 1);
////	cvShowImage("circles", img);
////	//cvShowImage("detected circles", &img);
////	cv::waitKey(0);
////
////
////	cvReleaseMemStorage(&storage);
////	cvReleaseImage(&imgGrayScale);
////	cvReleaseImage(&img);
////	
////	return true;
////}




