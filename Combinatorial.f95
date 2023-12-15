! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	Subroutine Combinatorial(iunit,np,npt,index,inum)
! 	Note that this subroutine is recursive
	integer*4 iunit,np,npt,index,inum
	integer*4 i1
	integer*4 in(20)
	common /combo/in
	
	
	do 10 i1 = inum,npt
	in(index) = i1
	if(np.eq.index)then
		write(iunit,*)(in(i),i=1,np)
		go to 10
		endif

	call Combinatorial(iunit,np,npt,index+1,i1+1)

10	continue

	return

	end

! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	Subroutine Combinatorial2(iunit,np,npt,index,inum)
! 	Note that this subroutine is recursive
	integer*4 iunit,np,npt,index,inum
	integer*4 i1
	integer*4 in(20)
	common /combo2/in
	
	
	do 10 i1 = inum,npt
	in(index) = i1
	if(np.eq.index)then
		write(iunit,*)np,(in(i),i=1,np)
		go to 10
		endif

	call Combinatorial2(iunit,np,npt,index+1,i1+1)

10	continue

	return

	end

