program leeweights    
    set output error
    version 16
    syntax varname, select(varname) TREATment(varname) [tight(varlist) suffix(name)]
    qui{
    if "`suffix'" != "" local suffix = "_`suffix'" 
    tempvar interaction
    tempvar s1_hat
    tempvar s0_hat
    tempvar tau_hat
    tempvar d_help
    tempvar d_hurt
    tempvar d_0
    tempvar tmp_rank_up
    tempvar tmp_rank_low
    tempvar rank_up
    tempvar rank_low
    tempvar count_n
    tempvar pct_utail
    tempvar pct_ltail
    tempvar ll
    tempvar ul
    // Generating survey response variable
    egen `interaction' = group(`tight')
    reg `select' i.(`interaction') if `treatment' == 1
    predict `s1_hat', xb
    reg `select' i.(`interaction') if `treatment' == 0
    predict `s0_hat', xb
    gen `tau_hat' = `s1_hat'-`s0_hat'
    gen `d_help'  = `tau_hat' > 0
    gen `d_hurt'  = `tau_hat' < 0
    gen `d_0'     = `tau_hat' == 0

    sort `tight' `treatment' `varlist'

    bys `tight' `treatment' (`varlist'): egen `tmp_rank_up' = rank(-`varlist') if !mi(`varlist'), unique
    bys `tight' `treatment' (`varlist'): egen `tmp_rank_low' = rank(`varlist') if !mi(`varlist'), unique
    bys `tight' `treatment' `varlist': egen `rank_up' = max(`tmp_rank_up') if !mi(`varlist')
    bys `tight' `treatment' `varlist': egen `rank_low' = max(`tmp_rank_low') if !mi(`varlist')
    bys `tight' `treatment' (`varlist'): egen `count_n' = count(`varlist') if !mi(`varlist')
    gen `pct_utail' = `rank_up'/`count_n' if !mi(`varlist')
    gen `pct_ltail' = `rank_low'/`count_n' if !mi(`varlist')

    *** Lee upper bound weights
    cap drop wgt_ub`suffix'
    gen wgt_ub`suffix' = 0
    replace wgt_ub`suffix' = 1 if `d_0' == 1
    ** Upper bound for X(help) - Trimming the lower tail of treatment
    replace wgt_ub`suffix' = 1 if `select' == 1 & `d_help' == 1 & `treatment' == 1 & `pct_ltail' > abs(`tau_hat') 
    * Reweighting for non-exact quantiles
    bys `tight' `treatment' (`varlist'): ///
        egen `ll' = ///
            max(0*(`pct_ltail'>=abs(`tau_hat')) + ///
            `pct_ltail'*(`pct_ltail'<abs(`tau_hat'))) ///
        if `select' == 1 & `d_help' == 1 & `treatment' == 1 
    bys `tight' `treatment' (`varlist'): ///
        egen `ul' = ///
            min(1*(`pct_ltail'<abs(`tau_hat')) + ///
            `pct_ltail'*(`pct_ltail'>=abs(`tau_hat'))) ///
        if `select' == 1 & `d_help' == 1 & `treatment' == 1 
    replace wgt_ub`suffix' = (`ul' - abs(`tau_hat')) / (`ul'-`ll') ///
        if `select' == 1 & `d_help' == 1 & `treatment' == 1 & `pct_ltail' > `ll' & `pct_ltail' <= `ul'  
    drop `ul' `ll'    

    ** Upper bound for X(help) - Full Control
    replace wgt_ub`suffix' = 1 if `select' == 1 & `d_help' == 1 & `treatment' == 0 

    ** Upper bound for X(hurt) - Full `treatment'
    replace wgt_ub`suffix' = 1 if `select' == 1 & `d_hurt' == 1 & `treatment' == 1

    ** Upper bound for X(hurt) - Trimming the upper tail of control
    replace wgt_ub`suffix' = 1 if `select' == 1 & `d_hurt' == 1 & `treatment' == 0 & `pct_utail' > abs(`tau_hat') 
    * Reweighting for non-exact quantiles
    bys `tight' `treatment' (`varlist'): ///
        egen `ll' = ///
            max(0*(`pct_utail'>=abs(`tau_hat')) + ///
            `pct_utail'*(`pct_utail'<abs(`tau_hat'))) ///
        if `select' == 1 & `d_hurt' == 1 & `treatment' == 0 
    bys `tight' `treatment' (`varlist'): ///
        egen `ul' = ///
            min(1*(`pct_utail'<abs(`tau_hat')) + ///
            `pct_utail'*(`pct_utail'>=abs(`tau_hat'))) ///
        if `select' == 1 & `d_hurt' == 1 & `treatment' == 0 
    replace wgt_ub`suffix' = (`ul' - abs(`tau_hat')) / (`ul'-`ll') ///
        if `select' == 1 & `d_hurt' == 1 & `treatment' == 0 & `pct_utail' > `ll' & `pct_utail' <= `ul'  
    drop `ul' `ll'    

    *** Lee lower bound weights
    cap drop wgt_lb`suffix'
    gen wgt_lb`suffix' = 0
    replace wgt_lb`suffix' = 1 if `d_0' == 1
    ** Lower bound for X(help) - Trimming the upper tail of treatment
    replace wgt_lb`suffix' = 1 if `select' == 1 & `d_help' == 1 & `treatment' == 1 & `pct_utail' > abs(`tau_hat') 
    * Reweighting when there are ties
    bys `tight' `treatment' (`varlist'): ///
        egen `ll' = ///
            max(0*(`pct_utail'>=abs(`tau_hat')) + ///
            `pct_utail'*(`pct_utail'<abs(`tau_hat'))) ///
        if `select' == 1 & `d_help' == 1 & `treatment' == 1 
    bys `tight' `treatment' (`varlist'): ///
        egen `ul' = ///
            min(1*(`pct_utail'<abs(`tau_hat')) + ///
            `pct_utail'*(`pct_utail'>=abs(`tau_hat'))) ///
        if `select' == 1 & `d_help' == 1 & `treatment' == 1 
    replace wgt_lb`suffix' = (`ul' - abs(`tau_hat')) / (`ul'-`ll') ///
        if `select' == 1 & `d_help' == 1 & `treatment' == 1 & `pct_utail' > `ll' & `pct_utail' <= `ul'  
    drop `ul' `ll'     

    ** Lower bound for X(help) - Full Control
    replace wgt_lb`suffix' = 1 if `select' == 1 & `d_help' == 1 & `treatment' == 0 

    ** Lower bound for X(hurt) - Full `treatment'
    replace wgt_lb`suffix' = 1 if `select' == 1 & `d_hurt' == 1 & `treatment' == 1

    ** Lower bound for X(hurt) - Trimming the lower tail of control
    replace wgt_lb`suffix' = 1 if `select' == 1 & `d_hurt' == 1 & `treatment' == 0 & `pct_ltail' > abs(`tau_hat') 
    * Reweighting for non-exact quantiles
    bys `tight' `treatment' (`varlist'): ///
        egen `ll' = ///
            max(0*(`pct_ltail'>=abs(`tau_hat')) + ///
            `pct_ltail'*(`pct_ltail'<abs(`tau_hat'))) ///
        if `select' == 1 & `d_hurt' == 1 & `treatment' == 0 
    bys `tight' `treatment' (`varlist'): ///
        egen `ul' = ///
            min(1*(`pct_ltail'<abs(`tau_hat')) + ///
            `pct_ltail'*(`pct_ltail'>=abs(`tau_hat'))) ///
        if `select' == 1 & `d_hurt' == 1 & `treatment' == 0 
    replace wgt_lb`suffix' = (`ul' - abs(`tau_hat')) / (`ul'-`ll') ///
        if `select' == 1 & `d_hurt' == 1 & `treatment' == 0 & `pct_ltail' > `ll' & `pct_ltail' <= `ul'  
    drop `ul' `ll'  
    drop `tmp_rank_up' `tmp_rank_low' `rank_up' `rank_low' `count_n' `pct_utail' `pct_ltail' 
    replace wgt_ub`suffix' = . if mi(`varlist')
    replace wgt_lb`suffix' = . if mi(`varlist')
end


