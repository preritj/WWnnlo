c    Author : Prerit Jaiswal
c    Subroutine to initialize the PDFs using LHAPDF subrotines 
c    When compiling, use -lLHApdf to include the PDF libraries 
c    PDF sets should be placed in usr/local/share/lhapdf/PDFsets/

******************************************************************
c      Subroutine pdfini()
c      Implicit None
c      Character nam*25
c      Integer mem             ! mem = member of PFset you want to use 
c      parameter(nam = "MRST2004nnlo.LHgrid")
c      parameter(mem = 0)      ! use mem = 0 for central PDF c

c      Call InitPDFsetbyname(nam)
c      Call InitPDF(mem)
      
c      Return
c      End

      Subroutine pdfini(Nloop)
      Implicit None
      Integer Nloop		! Loop order of PDF to be used
      Character nam*25
      Integer mem             ! mem = member of PFset you want to use 
       if (Nloop .eq. 0) then
    	nam = "MSTW2008lo68cl.LHgrid"
    	mem = 0
    	else if (Nloop .eq. 1) then
    	nam = "MSTW2008nlo68cl.LHgrid"
    	mem = 0
        else if (Nloop .gt. 1) then
         nam = "MSTW2008nnlo68cl.LHgrid"
         mem = 0
    	end if 
      Call InitPDFsetbyname(nam)
      Call InitPDF(mem)
      
      Return
      End

	
	Subroutine pdfevolve( x, Q, f)
	Implicit None
	double precision x, Q, f(-6:6)
	call evolvePDF(x, Q, f)
c	Print*, " FORTRAN SHIT  "
c	Print*, f(0), f(1), f(-1), f(4), f(-4)
	Return
	End
********************************************************************
