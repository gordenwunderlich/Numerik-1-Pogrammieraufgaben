module kinddef
    integer, parameter :: real_kind = 8
endmodule

module Tschebyscheff
    use kinddef
    contains
        pure function chebval(c, x) result(val)
            real(real_kind), dimension(:), intent(in) :: c, x
            real(real_kind), dimension(size(x)) :: val
            real(real_kind), dimension(size(x), 3) :: d
            d = 0
            do i = size(c), 1, -1
                d = reshape([c(i) + 2 * x * d(:, 2) - d(:, 3), d(:, 1), d(:, 2)],[size(x), 3])
            enddo
            val = d(:, 1)
        endfunction 
endmodule

program xycxcbdthfsfgb
    use kinddef
endprogram