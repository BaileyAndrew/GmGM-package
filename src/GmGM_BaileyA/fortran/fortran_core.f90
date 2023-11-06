module fortran_core
    implicit none

    contains
        subroutine extract_d_values(shape_, ndims, d_all, d_lefts, d_rights)
            ! shape represents the shape of a tensor,
            ! ndims is the number of dimensions, i.e. the
            ! length of shape.  for matrices it is 2.
            integer, intent(in) :: ndims
            integer, intent(in), dimension(ndims) :: shape_
        
            ! d_all: (3, 4, 6) -> 3 * 4 * 6 = 72
            ! d_lefts: (3, 4, 6) -> (1, 3, 12)
            ! d_rights: (3, 4, 6) -> (24, 6, 1)
            integer, intent(out) :: d_all
            integer, intent(out), dimension(ndims) :: d_lefts
            integer, intent(out), dimension(ndims) :: d_rights

            integer :: i
        
            do i = 1, ndims
                if (i == 1) then
                    d_lefts(i) = 1
                    d_rights(i) = product(shape_(2:ndims))
                else if (i == ndims) then
                    d_lefts(i) = product(shape_(1:i-1))
                    d_rights(i) = 1
                else
                    d_lefts(i) = product(shape_(1:i-1))
                    d_rights(i) = product(shape_(i+1:ndims))
                end if
            end do
        
            d_all = product(shape_)
        end subroutine

        ! TODO: combine into one
        subroutine project_inv_kron_sum_2(x, xs, y, ys, xproj, yproj)
            integer, intent(in) :: xs, ys
            real, intent(in), dimension(xs) :: x
            real, intent(in), dimension(ys) :: y
            real, intent(out), dimension(xs) :: xproj
            real, intent(out), dimension(ys) :: yproj
            
            real, parameter :: k_ratio = 1./2.
        
            real :: cur_val
            
            integer :: i, j
            
            do i=1,xs
                do j=1,ys
                    cur_val = 1 / (x(i)+y(j))
                    xproj(i) = xproj(i) + cur_val
                    yproj(j) = yproj(j) + cur_val
                end do
            end do
        
            ! Normalize
            xproj = xproj / ys
            yproj = yproj / xs
        
            ! Offset diagonal
            xproj = xproj - k_ratio * sum(xproj) / xs
            yproj = yproj - k_ratio * sum(yproj) / ys
            
        end subroutine
        
        subroutine project_inv_kron_sum_3(x, xs, y, ys, z, zs, xproj, yproj, zproj)
            integer, intent(in) :: xs, ys, zs
            real, intent(in), dimension(xs) :: x
            real, intent(in), dimension(ys) :: y
            real, intent(in), dimension(zs) :: z
            real, intent(out), dimension(xs) :: xproj
            real, intent(out), dimension(ys) :: yproj
            real, intent(out), dimension(zs) :: zproj
            
            real, parameter :: k_ratio = 2./3.
        
            real :: cur_val

            integer :: i, j, k
            
            do i=1,xs
                do j=1,ys
                    do k=1,zs
                        cur_val = 1 / (x(i)+y(j)+z(k))
                        xproj(i) = xproj(i) + cur_val
                        yproj(j) = yproj(j) + cur_val
                        zproj(k) = zproj(k) + cur_val
                    end do
                end do
            end do
        
            ! Normalize
            xproj = xproj / (ys*zs)
            yproj = yproj / (xs*zs)
            zproj = zproj / (xs*ys)
        
            ! Offset diagonal
            xproj = xproj - k_ratio * sum(xproj) / xs
            yproj = yproj - k_ratio * sum(yproj) / ys
            zproj = zproj - k_ratio * sum(zproj) / zs
        end subroutine
        
        subroutine project_inv_kron_sum_4(x, xs, y, ys, z, zs, w, ws, xproj, yproj, zproj, wproj)
            integer, intent(in) :: xs, ys, zs, ws
            real, intent(in), dimension(xs) :: x
            real, intent(in), dimension(ys) :: y
            real, intent(in), dimension(zs) :: z
            real, intent(in), dimension(ws) :: w
            real, intent(out), dimension(xs) :: xproj
            real, intent(out), dimension(ys) :: yproj
            real, intent(out), dimension(zs) :: zproj
            real, intent(out), dimension(ws) :: wproj
            
            real, parameter :: k_ratio = 3./4.
        
            real :: cur_val
            
            integer :: i, j, k, l

            do i=1,xs
                do j=1,ys
                    do k=1,zs
                        do l=1,ws
                            cur_val = 1 / (x(i)+y(j)+z(k)+w(l))
                            xproj(i) = xproj(i) + cur_val
                            yproj(j) = yproj(j) + cur_val
                            zproj(k) = zproj(k) + cur_val
                            wproj(l) = wproj(l) + cur_val
                        end do
                    end do
                end do
            end do
        
            ! Normalize
            xproj = xproj / (ys*zs*ws)
            yproj = yproj / (xs*zs*ws)
            zproj = zproj / (xs*ys*ws)
            wproj = wproj / (xs*ys*zs)
        
            ! Offset diagonal
            xproj = xproj - k_ratio * sum(xproj) / xs
            yproj = yproj - k_ratio * sum(yproj) / ys
            zproj = zproj - k_ratio * sum(zproj) / zs
            wproj = wproj - k_ratio * sum(wproj) / ws
        end subroutine
        
        subroutine project_inv_kron_sum_5(x, xs, y, ys, z, zs, w, ws, v, vs, xproj, yproj, zproj, wproj, vproj)
            integer, intent(in) :: xs, ys, zs, ws, vs
            real, intent(in), dimension(xs) :: x
            real, intent(in), dimension(ys) :: y
            real, intent(in), dimension(zs) :: z
            real, intent(in), dimension(ws) :: w
            real, intent(in), dimension(vs) :: v
            real, intent(out), dimension(xs) :: xproj
            real, intent(out), dimension(ys) :: yproj
            real, intent(out), dimension(zs) :: zproj
            real, intent(out), dimension(ws) :: wproj
            real, intent(out), dimension(vs) :: vproj
            
            real, parameter :: k_ratio = 4./5.
        
            real :: cur_val
            
            integer :: i, j, k, l, m

            do i=1,xs
                do j=1,ys
                    do k=1,zs
                        do l=1,ws
                            do m=1,vs
                                cur_val = 1 / (x(i)+y(j)+z(k)+w(l)+v(m))
                                xproj(i) = xproj(i) + cur_val
                                yproj(j) = yproj(j) + cur_val
                                zproj(k) = zproj(k) + cur_val
                                wproj(l) = wproj(l) + cur_val
                                vproj(m) = vproj(m) + cur_val
                            end do
                        end do
                    end do
                end do
            end do
        
            ! Normalize
            xproj = xproj / (ys*zs*ws*vs)
            yproj = yproj / (xs*zs*ws*vs)
            zproj = zproj / (xs*ys*ws*vs)
            wproj = wproj / (xs*ys*zs*vs)
            vproj = vproj / (xs*ys*zs*ws)
        
            ! Offset diagonal
            xproj = xproj - k_ratio * sum(xproj) / xs
            yproj = yproj - k_ratio * sum(yproj) / ys
            zproj = zproj - k_ratio * sum(zproj) / zs
            wproj = wproj - k_ratio * sum(wproj) / ws
            vproj = vproj - k_ratio * sum(vproj) / vs
        end subroutine
        
        subroutine project_inv_kron_sum_6(x, xs, y, ys, z, zs, w, ws, v, vs, u, us, xproj, yproj, zproj, wproj, vproj, uproj)
            integer, intent(in) :: xs, ys, zs, ws, vs, us
            real, intent(in), dimension(xs) :: x
            real, intent(in), dimension(ys) :: y
            real, intent(in), dimension(zs) :: z
            real, intent(in), dimension(ws) :: w
            real, intent(in), dimension(vs) :: v
            real, intent(in), dimension(us) :: u
            real, intent(out), dimension(xs) :: xproj
            real, intent(out), dimension(ys) :: yproj
            real, intent(out), dimension(zs) :: zproj
            real, intent(out), dimension(ws) :: wproj
            real, intent(out), dimension(vs) :: vproj
            real, intent(out), dimension(us) :: uproj
            
            real, parameter :: k_ratio = 5./6.
        
            real :: cur_val

            integer :: i, j, k, l, m, n
            
            do i=1,xs
                do j=1,ys
                    do k=1,zs
                        do l=1,ws
                            do m=1,vs
                                do n=1,us
                                    cur_val = 1 / (x(i)+y(j)+z(k)+w(l)+v(m)+u(n))
                                    xproj(i) = xproj(i) + cur_val
                                    yproj(j) = yproj(j) + cur_val
                                    zproj(k) = zproj(k) + cur_val
                                    wproj(l) = wproj(l) + cur_val
                                    vproj(m) = vproj(m) + cur_val
                                    uproj(n) = uproj(n) + cur_val
                                end do
                            end do
                        end do
                    end do
                end do
            end do
        
            ! Normalize
            xproj = xproj / (ys*zs*ws*vs*us)
            yproj = yproj / (xs*zs*ws*vs*us)
            zproj = zproj / (xs*ys*ws*vs*us)
            wproj = wproj / (xs*ys*zs*vs*us)
            vproj = vproj / (xs*ys*zs*ws*us)
            uproj = uproj / (xs*ys*zs*ws*vs)
        
            ! Offset diagonal
            xproj = xproj - k_ratio * sum(xproj) / xs
            yproj = yproj - k_ratio * sum(yproj) / ys
            zproj = zproj - k_ratio * sum(zproj) / zs
            wproj = wproj - k_ratio * sum(wproj) / ws
            vproj = vproj - k_ratio * sum(vproj) / vs
            uproj = uproj - k_ratio * sum(uproj) / us
        end subroutine

        subroutine sum_log_sum_2(x, xs, y, ys, sl)
            integer, intent(in) :: xs, ys
            real, dimension(xs), intent(in) :: x
            real, dimension(ys), intent(in) :: y
            real, intent(out) :: sl
            
            integer :: idx = 1
            real :: intermediate = 1
            integer :: simplify_size = 20
            real :: sm = 0
        
            integer :: i, j

            sl = 0
            
            do i=1,xs
                do j=1,ys
                    sm = x(i)+y(j)
                    intermediate = intermediate * sm
                    if (idx == simplify_size) then
                        sl = sl + log(intermediate)
                        intermediate = 1.0
                        idx = 1
                    end if
                    idx = idx + 1
                end do
            end do
        
            sl = sl + log(intermediate)
        end subroutine
        
        subroutine sum_log_sum_3(x, xs, y, ys, z, zs, sl)
            integer, intent(in) :: xs, ys, zs
            real, dimension(xs), intent(in) :: x
            real, dimension(ys), intent(in) :: y
            real, dimension(zs), intent(in) :: z
            real, intent(out) :: sl
            
            integer :: idx = 1
            real :: intermediate = 1
            integer :: simplify_size = 20

            integer :: i, j, k
            
            do i=1,xs
                do j=1,ys
                    do k=1,zs
                        intermediate = intermediate * (x(i)+y(j)+z(k))
                        if (idx == simplify_size) then
                            sl = sl + log(intermediate)
                            intermediate = 1
                            idx = 1
                        end if
                        idx = idx + 1
                    end do
                end do
            end do
            return
        end subroutine
        
        subroutine sum_log_sum_4(x, xs, y, ys, z, zs, w, ws, sl)
            integer, intent(in) :: xs, ys, zs, ws
            real, dimension(xs), intent(in) :: x
            real, dimension(ys), intent(in) :: y
            real, dimension(zs), intent(in) :: z
            real, dimension(ws), intent(in) :: w
            real, intent(out) :: sl
            
            integer :: idx = 1
            real :: intermediate = 1
            integer :: simplify_size = 20

            integer :: i, j, k, l
            
            do i=1,xs
                do j=1,ys
                    do k=1,zs
                        do l=1,ws
                            intermediate = intermediate * (x(i)+y(j)+z(k)+w(l))
                            if (idx == simplify_size) then
                                sl = sl + log(intermediate)
                                intermediate = 1
                                idx = 1
                            end if
                            idx = idx + 1
                        end do
                    end do
                end do
            end do
            sl = sl + log(intermediate)
        end subroutine
        
        subroutine sum_log_sum_5(x, xs, y, ys, z, zs, w, ws, v, vs, sl)
            integer, intent(in) :: xs, ys, zs, ws, vs
            real, dimension(xs), intent(in) :: x
            real, dimension(ys), intent(in) :: y
            real, dimension(zs), intent(in) :: z
            real, dimension(ws), intent(in) :: w
            real, dimension(vs), intent(in) :: v
            real, intent(out) :: sl
            
            integer :: idx = 1
            real :: intermediate = 1
            integer :: simplify_size = 20

            integer :: i, j, k, l, m
            
            do i=1,xs
                do j=1,ys
                    do k=1,zs
                        do l=1,ws
                            do m=1,vs
                                intermediate = intermediate * (x(i)+y(j)+z(k)+w(l)+v(m))
                                if (idx == simplify_size) then
                                    sl = sl + log(intermediate)
                                    intermediate = 1
                                    idx = 1
                                end if
                                idx = idx + 1
                            end do
                        end do
                    end do
                end do
            end do
            sl = sl + log(intermediate)
        end subroutine
        
        subroutine sum_log_sum_6(x, xs, y, ys, z, zs, w, ws, v, vs, u, us, sl)
            integer, intent(in) :: xs, ys, zs, ws, vs, us
            real, dimension(xs), intent(in) :: x
            real, dimension(ys), intent(in) :: y
            real, dimension(zs), intent(in) :: z
            real, dimension(ws), intent(in) :: w
            real, dimension(vs), intent(in) :: v
            real, dimension(us), intent(in) :: u
            real, intent(out) :: sl
            
            integer :: idx = 1
            real :: intermediate = 1
            integer :: simplify_size = 20

            integer :: i, j, k, l, m, n
            
            do i=1,xs
                do j=1,ys
                    do k=1,zs
                        do l=1,ws
                            do m=1,vs
                                do n=1,us
                                    intermediate = intermediate * (x(i)+y(j)+z(k)+w(l)+v(m)+u(n))
                                    if (idx == simplify_size) then
                                        sl = sl + log(intermediate)
                                        intermediate = 1
                                        idx = 1
                                    end if
                                    idx = idx + 1
                                end do
                            end do
                        end do
                    end do
                end do
            end do
            sl = sl + log(intermediate)
        end subroutine

end module fortran_core