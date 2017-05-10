module kinddef
    integer, parameter :: real_kind = 8
endmodule

module mittelpunkt_regel
        use kinddef
    interface
        elemental function func(x) result(res)
            use kinddef
            real(real_kind), intent(in) :: x
            real(real_kind) :: res
        endfunction
    endinterface
    
    contains
        function mittelp(f, a, b, n) result(approx)
            real(real_kind), intent(in) :: a, b
            real(real_kind) :: approx
            integer, intent(in) :: n
            procedure(func) :: f
            real(real_kind) :: h
            h = (b-a)/n
            approx = sum( f((/((i+0.5_8)*h+a,i = 0,N-1)/)) )*h
        endfunction
endmodule

program xycxcbdthfsfgb
    use mittelpunkt_regel, only : mittelp
    use kinddef    
    real(real_kind) :: x
    real(real_kind), intrinsic :: dcos
    x = mittelp(dcos, 0._8, 1._8, 10)
    print *, x, "Fehler: ", x-sin(1._8)
    x = mittelp(dcos, 0._8, 1._8, 100)
    print *, x, "Fehler: ", x-sin(1._8)
    x = mittelp(dcos, 0._8, 1._8, 1000)
    print *, x, "Fehler: ", x-sin(1._8)
endprogram