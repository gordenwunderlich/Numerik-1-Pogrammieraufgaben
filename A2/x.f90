module kinddef
    integer, parameter :: real_kind = 8
endmodule

module Quadratur
        use kinddef
    interface
        elemental function func(x) result(res)
            use kinddef
            real(real_kind), intent(in) :: x
            real(real_kind) :: res
        endfunction
    endinterface
    
    contains
        pure function numint(f, al, be, n, b, c) result(intapprox)
            real(real_kind), intent(in) :: al, be
            real(real_kind) :: intapprox
            real(real_kind), dimension(:), intent(in) :: b, c
            integer, intent(in) :: n
            procedure(func) :: f
            real(real_kind) :: h
            h = (be - al) / N
            intapprox = sum(f( (/((i+c)*h+al, i=0,N-1)/) )* reshape(spread(b,2,N),(/30/)) * h)
        endfunction
endmodule

program xycxcbdthfsfgb
    use Quadratur, only : numint
    use omp_lib, only : omp_get_wtime
    use kinddef    
    character(1), dimension(6) :: number
    real(real_kind), dimension(7,6) :: x
    real(real_kind), dimension(7, 6) :: timestart, timeend
    real(real_kind), parameter :: exakt = 1._8 - cos(1._8)
    integer, dimension(7), parameter :: N = (/ 10, 20, 50, 100, 200, 500, 1000/)
    !real(real_kind), intrinsic :: dcos
    do i=1,7
        timestart(i, 1) = omp_get_wtime()
        x(i,1) = numint(func, 0._real_kind, 1._real_kind, N(i), (/ 1._real_kind /), (/ 0._real_kind /)) !rechteckregel
        timeend(i, 1) = omp_get_wtime()
    enddo
    do i=1,7
        timestart(i, 2) = omp_get_wtime()
        x(i,2) = numint(func, 0._real_kind, 1._real_kind, N(i), (/ 1._real_kind /), (/ 0.5_real_kind /)) !mittelpunktregel
        timeend(i, 2) = omp_get_wtime()
    enddo
    do i=1,7
        timestart(i, 3) = omp_get_wtime()
        x(i,3) = numint(func, 0._real_kind, 1._real_kind, N(i), (/ 0.5_real_kind, 0.5_real_kind /), (/ 0._real_kind, 1._real_kind /)) !trapezregel
        timeend(i, 3) = omp_get_wtime()
    enddo
    do i=1,7
        timestart(i, 4) = omp_get_wtime()
        x(i,4) = numint(func, 0._real_kind, 1._real_kind, N(i), (/ 1._real_kind/6._real_kind, 4._real_kind/6._real_kind, 1._real_kind/6._real_kind /), (/ 0._real_kind, 0.5_real_kind, 1._real_kind /)) !simpsonregel
        timeend(i, 4) = omp_get_wtime()
    enddo
    do i=1,7
        timestart(i, 5) = omp_get_wtime()
        x(i,5) = numint(func, 0._real_kind, 1._real_kind, N(i), (/ 0.75_real_kind, 0.25_real_kind /), (/ 1._real_kind/3._real_kind, 1._real_kind /)) !verfahren 1
        timeend(i, 5) = omp_get_wtime()
    enddo
    do i=1,7
        timestart(i, 6) = omp_get_wtime()
        x(i,6) = numint(func, 0._real_kind, 1._real_kind, N(i), (/ 0.5_real_kind, 0.5_real_kind /), (/ 0.5_real_kind - (sqrt(3._real_kind)/6._real_kind ), 0.5_real_kind + (sqrt(3._real_kind)/6._real_kind ) /)) !verfahren 2
        timeend(i, 6) = omp_get_wtime()
    enddo
    timeend = timeend
    print "(A)", "     Rechteckregel                           Mittelpunktregel                       Trapezregel"
    print "(7(3(5x ,A, i4, a, g0), /, 3(5x, a, g0, 8x), /))", ("N = ", N(i), "; F: ", abs(x(i,1) - exakt), "N = ", N(i), "; F: ", abs(x(i,2) - exakt), "N = ", N(i), "; F: ", abs(x(i,3) - exakt), ("Time: ", timeend(i, j) - timestart(i, j), j = 1,3), i = 1, 7)
    print "(A)", "     Simpsonregel                           Verfahren 1                             Verfahren 2"
    print "(7(3(5x ,A, i4, a, g0), /, 3(5x, a, g0, 7x), /))", ("N = ", N(i), "; F: ", abs(x(i,4) - exakt), "N = ", N(i), "; F: ", abs(x(i,5) - exakt), "N = ", N(i), "; F: ", abs(x(i,6) - exakt), ("Time: ", timeend(i, j) - timestart(i, j), j = 4,6), i = 1, 7)
    write(number, "(i1)") (j, j = 1,6)
    do concurrent( j = 1: 6)
        open(newunit = iunit, file="dataN"// number(j) //".dat", action="write")
        write(iunit, "(*(g0, x, g0, /))") (N(i),abs(x(i,j)- exakt), i = 1,7) 
        close(iunit)
    enddo
    
    do concurrent( j = 1: 6)
        open(newunit = iunit, file="dataT"// number(j) //".dat", action="write")
        write(iunit, "(*(g0, x, g0, /))") ((timeend(i, j) - timestart(i, j)),abs(x(i,j)- exakt), i = 1,7) 
        close(iunit)
    enddo
    
        call System("gnuplot /Users/Gorden/Desktop/Numerik/A2/plotscript.txt")
    contains
        elemental function func(x) result(res)
            use kinddef
            use omp_lib
            real(real_kind), intent(in) :: x
            real(real_kind) :: res
            real(8) :: start, counter
            start = omp_get_wtime()
            counter = 0.
            do while(counter .lt. start + 0.000025)
                counter = omp_get_wtime()
            enddo
            res = sin(x)
        endfunction
endprogram