!    Author : Prerit Jaiswal
!    Subroutine to initialize the PDFs using LHAPDF subroutines 
!    When compiling, use -lLHApdf to include the PDF libraries 

!******************************************************************
      Subroutine pdfini(nam, mem)
      Implicit None
      Character*35 nam        ! nam = name of the PDF set 
      Integer mem             ! mem = member of PDF set 

      Call InitPDFsetbyname(nam)
!      Call InitPDFsetbyname('MSTW2008nnlo90cl_nf4.LHgrid')
      Call InitPDF(mem)
!      Call InitPDF(0)
     
      Return
      End


!******************************************************************

	
	Subroutine pdfevolve( x, Q, f)
	Implicit None
	double precision x, Q, f(-6:6)
	call evolvePDF(x, Q, f)
	Return
	End
!********************************************************************
