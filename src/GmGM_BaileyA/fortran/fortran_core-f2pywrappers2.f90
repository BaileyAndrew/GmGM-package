!     -*- f90 -*-
!     This file is autogenerated with f2py (version:1.22.4)
!     It contains Fortran 90 wrappers to fortran functions.

      
      subroutine f2pyinitfortran_core(f2pysetupfunc)
      use fortran_core, only : extract_d_values
      use fortran_core, only : project_inv_kron_sum_2
      use fortran_core, only : project_inv_kron_sum_3
      use fortran_core, only : project_inv_kron_sum_4
      use fortran_core, only : project_inv_kron_sum_5
      use fortran_core, only : project_inv_kron_sum_6
      use fortran_core, only : sum_log_sum_2
      use fortran_core, only : sum_log_sum_3
      use fortran_core, only : sum_log_sum_4
      use fortran_core, only : sum_log_sum_5
      use fortran_core, only : sum_log_sum_6
      external f2pysetupfunc
      call f2pysetupfunc(extract_d_values,project_inv_kron_sum_2,project&
     &_inv_kron_sum_3,project_inv_kron_sum_4,project_inv_kron_sum_5,proj&
     &ect_inv_kron_sum_6,sum_log_sum_2,sum_log_sum_3,sum_log_sum_4,sum_l&
     &og_sum_5,sum_log_sum_6)
      end subroutine f2pyinitfortran_core

