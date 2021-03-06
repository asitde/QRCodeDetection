
// QRCodeDetectionDlg.cpp : implementation file
//

#include "stdafx.h"
#include "QRCodeDetection.h"
#include "QRCodeDetectionDlg.h"
#include "afxdialogex.h"
#include "QRCodeDetector.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#endif


// CAboutDlg dialog used for App About

class CAboutDlg : public CDialogEx
{
public:
	CAboutDlg();

// Dialog Data
#ifdef AFX_DESIGN_TIME
	enum { IDD = IDD_ABOUTBOX };
#endif

	protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support

// Implementation
protected:
	DECLARE_MESSAGE_MAP()
};

CAboutDlg::CAboutDlg() : CDialogEx(IDD_ABOUTBOX)
{
}

void CAboutDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialogEx::DoDataExchange(pDX);
}

BEGIN_MESSAGE_MAP(CAboutDlg, CDialogEx)
END_MESSAGE_MAP()


// CQRCodeDetectionDlg dialog



CQRCodeDetectionDlg::CQRCodeDetectionDlg(CWnd* pParent /*=NULL*/)
	: CDialog(IDD_QRCODEDETECTION_DIALOG, pParent)
{
	m_hIcon = AfxGetApp()->LoadIcon(IDR_MAINFRAME);
}

void CQRCodeDetectionDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
}

BEGIN_MESSAGE_MAP(CQRCodeDetectionDlg, CDialog)
	ON_WM_SYSCOMMAND()
	ON_WM_PAINT()
	ON_WM_QUERYDRAGICON()
	ON_BN_CLICKED(IDC_BTN_LOAD_IMAGE, &CQRCodeDetectionDlg::OnBnClickedBtnLoadImage)
	//ON_MESSAGE(QRCodeDetector::UWM_DISP_DISTANCE, OnDispDistance)
	//ON_MESSAGE(1001, OnDispDistance)
END_MESSAGE_MAP()


// CQRCodeDetectionDlg message handlers

BOOL CQRCodeDetectionDlg::OnInitDialog()
{
	CDialog::OnInitDialog();

	// Add "About..." menu item to system menu.

	// IDM_ABOUTBOX must be in the system command range.
	ASSERT((IDM_ABOUTBOX & 0xFFF0) == IDM_ABOUTBOX);
	ASSERT(IDM_ABOUTBOX < 0xF000);

	CMenu* pSysMenu = GetSystemMenu(FALSE);
	if (pSysMenu != NULL)
	{
		BOOL bNameValid;
		CString strAboutMenu;
		bNameValid = strAboutMenu.LoadString(IDS_ABOUTBOX);
		ASSERT(bNameValid);
		if (!strAboutMenu.IsEmpty())
		{
			pSysMenu->AppendMenu(MF_SEPARATOR);
			pSysMenu->AppendMenu(MF_STRING, IDM_ABOUTBOX, strAboutMenu);
		}
	}

	// Set the icon for this dialog.  The framework does this automatically
	//  when the application's main window is not a dialog
	SetIcon(m_hIcon, TRUE);			// Set big icon
	SetIcon(m_hIcon, FALSE);		// Set small icon

	// TODO: Add extra initialization here
	SetDlgItemText(IDC_EDIT_MIN_RADIUS, _T("650"));
	SetDlgItemText(IDC_EDIT_MAX_RADIUS, _T("765"));
	SetDlgItemText(IDC_EDIT_PARAM1, _T("200"));
	SetDlgItemText(IDC_EDIT_PARAM2, _T("20"));


	return TRUE;  // return TRUE  unless you set the focus to a control
}

void CQRCodeDetectionDlg::OnSysCommand(UINT nID, LPARAM lParam)
{
	if ((nID & 0xFFF0) == IDM_ABOUTBOX)
	{
		CAboutDlg dlgAbout;
		dlgAbout.DoModal();
	}
	else
	{
		CDialog::OnSysCommand(nID, lParam);
	}
}

// If you add a minimize button to your dialog, you will need the code below
//  to draw the icon.  For MFC applications using the document/view model,
//  this is automatically done for you by the framework.

void CQRCodeDetectionDlg::OnPaint()
{
	if (IsIconic())
	{
		CPaintDC dc(this); // device context for painting

		SendMessage(WM_ICONERASEBKGND, reinterpret_cast<WPARAM>(dc.GetSafeHdc()), 0);

		// Center icon in client rectangle
		int cxIcon = GetSystemMetrics(SM_CXICON);
		int cyIcon = GetSystemMetrics(SM_CYICON);
		CRect rect;
		GetClientRect(&rect);
		int x = (rect.Width() - cxIcon + 1) / 2;
		int y = (rect.Height() - cyIcon + 1) / 2;

		// Draw the icon
		dc.DrawIcon(x, y, m_hIcon);
	}
	else
	{
		CDialog::OnPaint();
	}
}

// The system calls this function to obtain the cursor to display while the user drags
//  the minimized window.
HCURSOR CQRCodeDetectionDlg::OnQueryDragIcon()
{
	return static_cast<HCURSOR>(m_hIcon);
}



void CQRCodeDetectionDlg::OnBnClickedBtnLoadImage()
{

	CString csRoot{ _T("C:\\Images") };
	QRCodeDetector myDetector;

	CString csMinRadius{ _T("") };
	GetDlgItemText(IDC_EDIT_MIN_RADIUS, csMinRadius);
	std::string strMinRadius{ CT2CA{csMinRadius}.m_psz };
	int iMinRadius = std::stoi(strMinRadius);

	CString csMaxRadius{ _T("") };
	GetDlgItemText(IDC_EDIT_MAX_RADIUS, csMaxRadius);
	std::string strMaxRadius{ CT2CA{ csMaxRadius }.m_psz };
	int iMaxRadius = std::stoi(strMaxRadius);
	myDetector.SetRadii( iMinRadius, iMaxRadius);

	CString csParam1{ _T("") };
	GetDlgItemText(IDC_EDIT_PARAM1, csParam1);
	std::string strParam1{ CT2CA{ csParam1 }.m_psz };
	double dParam1 = std::stof(strParam1);

	CString csParam2{ _T("") };
	GetDlgItemText(IDC_EDIT_PARAM2, csParam2);
	std::string strParam2{ CT2CA{ csParam2 }.m_psz };
	double dParam2 = std::stof(strParam2);
	myDetector.SetParams(dParam1, dParam2);

	GetDlgItem(IDC_EDIT_MIN_RADIUS)->EnableWindow(FALSE);
	GetDlgItem(IDC_EDIT_MAX_RADIUS)->EnableWindow(FALSE);
	GetDlgItem(IDC_EDIT_PARAM1)->EnableWindow(FALSE);
	GetDlgItem(IDC_EDIT_PARAM2)->EnableWindow(FALSE);


	myDetector.Run(csRoot);

	GetDlgItem(IDC_EDIT_MIN_RADIUS)->EnableWindow(TRUE);
	GetDlgItem(IDC_EDIT_MAX_RADIUS)->EnableWindow(TRUE);
	GetDlgItem(IDC_EDIT_PARAM1)->EnableWindow(TRUE);
	GetDlgItem(IDC_EDIT_PARAM2)->EnableWindow(TRUE);

}
