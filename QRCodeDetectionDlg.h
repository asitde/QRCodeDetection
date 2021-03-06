
// QRCodeDetectionDlg.h : header file
//

#pragma once
#include "afxwin.h"


// CQRCodeDetectionDlg dialog
class CQRCodeDetectionDlg : public CDialog
{
// Construction
public:
	CQRCodeDetectionDlg(CWnd* pParent = NULL);	// standard constructor

// Dialog Data
#ifdef AFX_DESIGN_TIME
	enum { IDD = IDD_QRCODEDETECTION_DIALOG };
#endif

	protected:
	virtual void DoDataExchange(CDataExchange* pDX);	// DDX/DDV support


// Implementation
protected:
	HICON m_hIcon;

	// Generated message map functions
	virtual BOOL OnInitDialog();
	afx_msg void OnSysCommand(UINT nID, LPARAM lParam);
	afx_msg void OnPaint();
	afx_msg HCURSOR OnQueryDragIcon();
	DECLARE_MESSAGE_MAP()
public:
	afx_msg void OnBnClickedBtnLoadImage();
	////afx_msg LRESULT OnDispDistance(WPARAM, LPARAM);
private:
	CEdit edL_;
	CEdit edR_;
	CEdit edT_;
	CEdit edB_;
};
