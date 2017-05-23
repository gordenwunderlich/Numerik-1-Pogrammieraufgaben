module kinddef
    integer, parameter :: real_kind = 8
endmodule

module Fassregel
        use kinddef
        real(16), parameter :: PI = 4*atan(1._16)
        type kepler_result
            real(real_kind) :: vol
            real(real_kind), dimension(:,:), allocatable :: points
        endtype
        contains
            pure function kepler(phi, l) result(res)
                real(real_kind), dimension(:), intent(in) :: phi, l
                real(real_kind), dimension(:), allocatable :: B, H
                real(real_kind) ::  A
                type(kepler_result) :: res
                A = l(size(l)/2+1)
                B = (cos(phi) * l - A/2)
                H = sin(phi) * l
                res%vol = sum((B(:size(B) - 1) + B(2:)) / 2._real_kind * (h(:size(h) - 1) - h(2:))) ** 2 * PI / 2
                res%points = reshape([B,H],[2,size(b)], order = [2,1])
            endfunction
endmodule



module Fassdaten
    use kinddef
    real(real_kind), dimension(*), parameter :: phi1 = [&
  0.4636476090008061 ,& 
  0.3708912888126626 ,& 
  0.2709468503384208 ,& 
  0.1651486774146261 ,& 
  0.0554985052457156 ,& 
  0. ,&
  -0.0554985052457156 ,& 
  -0.1651486774146261 ,& 
  -0.2709468503384208 ,& 
  -0.3708912888126626 ,& 
  -0.4636476090008061 & 
], &
l1 = [&
  2.2360679774997898 ,& 
  2.1459119906475519 ,& 
  2.0757268546966006 ,& 
  2.0275875100994063 ,& 
  2.0030840419244385 ,& 
  2. ,&
  2.0030840419244385 ,& 
  2.0275875100994063 ,& 
  2.0757268546966006 ,& 
  2.1459119906475519 ,& 
  2.2360679774997898 & 
], &
phi2 = [&
  0.5070985043923368 ,& 
  0.3924562028632983 ,& 
  0.2791407016999940 ,& 
  0.1669701690249294 ,& 
  0.0555669655728478 ,& 
  0. ,&
  -0.0555669655728478 ,& 
  -0.1669701690249294 ,& 
  -0.2791407016999940 ,& 
  -0.3924562028632983 ,& 
  -0.5070985043923368 & 
], &
l2 = [&
  2.0591260281974000 ,& 
  2.0336237771080183 ,& 
  2.0163181271363468 ,& 
  2.0056709723637915 ,& 
  2.0006187124072592 ,& 
  2. ,&
  2.0006187124072592 ,& 
  2.0056709723637915 ,& 
  2.0163181271363468 ,& 
  2.0336237771080183 ,& 
  2.0591260281974000 & 
]
real(real_kind), parameter :: ext1 = 6.2831853071795862, ext2 = 5.4956927486797440


endmodule

program xycxcbdthfsfgb
    use Fassregel, only : kepler, kepler_result
    use omp_lib, only : omp_get_wtime
    use kinddef
    use Fassdaten
    type(kepler_result) :: res1, res2
    res1 = kepler(phi1, l1)
    print *, res1%vol
    print *, ext1
    res2 = kepler(phi2, l2)
    print *, res2%vol
    print *, ext2
    open(newunit = i, file = "")
    write "(*(2f10.5, /))", i, res1%points
endprogram