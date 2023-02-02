!> \file mo_spem.f90

!> \brief Routines for bias insensitive comparison of spatial patterns.

!> \author Moctar Dembélé
!> \date Dec 2018

MODULE mo_spem

  USE mo_kind, ONLY : i4, sp, dp

  IMPLICIT NONE

  PUBLIC :: spem       ! number of neighboring dominating values

  ! ------------------------------------------------------------------

  !     NAME
  !         SPEM - Spatial Pattern Efficiency Metric

  !     PURPOSE
  !>         \brief    Calculates the number of neighboring dominating values, a measure for spatial dissimilarity.
  !>         \details
  !>             SPEM = 1-sqrt((1-alpha)^2 + (1-beta)^2 + (1-gamma)^2)
  !>             alpha = corr(obs_vec,sim_vec,'rows','pairwise','type','Spearman');%Spearman correlation coefficient
  !>             beta = (nanstd(sim_vec)*nanmean(obs_vec))/(nanstd(obs_vec)*nanmean(sim_vec));%ratio of coefficient of variation
  !>             gamma = 1-sqrt(nanmean((zs_sim-zs_obs).^2));%1-RMSE of z-scores

  !     CALLING SEQUENCE
  !         out = SPEM(vec1, vec2, mask=mask)

  INTERFACE spem
    MODULE PROCEDURE spem_dp
  END INTERFACE spem

  ! ------------------------------------------------------------------

CONTAINS

  FUNCTION spem_dp(obs, sim, mask)
  
    use mo_errormeasures, only : rmse
    use mo_message, only : message
    use mo_common_constants, only : eps_dp, nodata_dp
    use mo_moment, only : mean, stddev, spearman_rho
    use mo_standard_score, only : standard_score

    IMPLICIT NONE

    real(dp), dimension(:), intent(IN) :: obs, sim
    logical, dimension(:), optional, intent(IN) :: mask
    logical, dimension(size(obs)) :: maske
    real(dp) :: n, spem_dp
    real(dp) :: alpha, beta, gamma
    real(dp) :: mean_obs, stddev_obs, cv_obs
    real(dp) :: mean_sim, stddev_sim, cv_sim
    real(dp), dimension(:), allocatable :: obs_final, sim_final
    real(dp), dimension(:), allocatable :: z_score_obs, z_score_sim
  
    !print *, 'In mo_spem: size(obs) = ', size(obs)
    !print *, 'In mo_spem: size(sim) = ', size(sim)
    if (size(obs) .ne. size(sim)) stop 'Error spem_dp: size(obs) .ne. size(sim)'
    if (present(mask)) then
      if (size(mask) .ne. size(obs)) stop 'Error spem_dp: size(mask) .ne. size(obs)'
      maske = mask
      n = real(count(maske), dp)
    else
      maske(:) = .true.
      n = real(size(obs), dp)
    end if
    if (n .le. (1.0_dp + tiny(1.0_dp))) stop 'spem_dp: n must be at least 2'
	
    allocate(z_score_obs(n))
    allocate(z_score_sim(n))
    z_score_obs(:) = nodata_dp
    z_score_sim(:) = nodata_dp
	
    !define mask for missing data in observations (z-score calculation does not tolerate nodata)
    obs_final = obs
    sim_final = sim
    where(abs(obs_final - nodata_dp) .lt. eps_dp) sim_final = nodata_dp
    where(abs(sim_final - nodata_dp) .lt. eps_dp) obs_final = nodata_dp
    where(abs(obs_final - nodata_dp) .lt. eps_dp) maske = .FALSE.
	
	!Calculate alpha
    alpha = spearman_rho(obs_final, sim_final, maske)
	
	!Calculate beta (ratio of coefficient of variation)
    mean_obs = mean(obs_final, maske)
    mean_sim = mean(sim_final, maske)
    stddev_obs = stddev(obs_final, maske)
    stddev_sim = stddev(sim_final, maske)
	!To prevent errors (i.e. denominator cannot be 0) 
    if (stddev_obs == 0_dp .OR. stddev_sim == 0_dp .OR. mean_obs == 0_dp .OR. mean_sim == 0_dp) then
      beta = 0_dp
      call message('WARNING: beta term in SPEM was forced to be 0')
      call message('i.e. at least stddev_obs, mean_obs, stddev_sim or mean_sim is 0)')
    else
      cv_obs = stddev_obs / mean_obs
      cv_sim = stddev_sim / mean_sim
      beta = cv_sim / cv_obs
    end if
	
	!Calculate gamma (RMSE (z-score_sim, z-score_obs))
    z_score_sim = standard_score(sim_final, maske)
    z_score_obs = standard_score(obs_final, maske)
    gamma = 1 - rmse(z_score_sim, z_score_obs, maske)
	
	!Calculate SPEM (Spatial Pattern Efficiency Metric)
    spem_dp = 1 - sqrt((1- alpha)**2 + (1 - beta)**2 + (1 - gamma)**2)
	
  END FUNCTION spem_dp

  ! ----------------------------------------------------------------------------------------------------------------

END MODULE mo_spem
