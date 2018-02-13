#pragma once

#include <string>
#include "highgui.h"
#include "opencv2/imgcodecs.hpp"
#include "opencv2/imgproc/imgproc.hpp"
#include <opencv2/opencv.hpp>
#include <fstream>
#include <cmath>

class QRCodeDetector {

public:

	struct QRDistance
	{
		int iLeft_;
		int iRight_;
		int iTop_;
		int iBottom_;
	};

	CONST  UINT UWM_DISP_DISTANCE = 1001;

	// Constructor
	QRCodeDetector();

	~QRCodeDetector();

	QRCodeDetector(const QRCodeDetector& rhs) = delete;

	QRCodeDetector& operator = (const QRCodeDetector& rhs) = delete;

	bool Run(const CString& rcstrImageDirectory);

	void SetRadii(int iMinRadius, int iMaxRadius);

	void SetParams(double dParam1, double dParam2);

private:

	bool Process(const std::string& rstrImageFile);
	bool ProcessEx(const std::string& rstrImageFile);

	void ApplyGaussianBlur(const cv::Mat& rOriginalImage, cv::Mat& rResultantImage);
	void ApplyMedianBlur(const cv::Mat& rOriginalImage, cv::Mat& rResultantImage);
	void DetectCircles(const cv::Mat& rSrcImage, std::vector<cv::Vec3f>& rVecCircles);
	void ConvertToGrayScale(const cv::Mat& rOriginalImage, cv::Mat& rResultantImage);
	void ApplyThresold(const cv::Mat& rOriginalImage, cv::Mat& rResultantImage);
	void ApplyAdaptiveThresold(const cv::Mat& rOriginalImage, cv::Mat& rResultantImage);
	void DrawDetectedCircles(cv::Mat& rOriginalImage, const std::vector<cv::Vec3f>& rVecCircles);
	void ShowImage(const cv::Mat& rSrcImage, const std::string& rstrWindowName);
	void SearchContours(const cv::Mat& rSrcImage, std::vector<std::vector<cv::Point>>& rVecContours);

public:

	static const double PIXEL_TO_MM;

private:

	std::string		m_strImagePath;
	std::string		m_strRootDirectory;
	std::string		m_strOutputDirectory;
	std::ofstream	m_ofs;
	double			m_dParam1;
	double			m_dParam2;
	int				m_iMinRadius;
	int				m_iMaxRadius;
};