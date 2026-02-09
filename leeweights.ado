program leeweights    
    set output error
    version 16
    syntax varname [if/] [pweight aweight/], select(varname) TREATment(varname) [tight(varlist) suffix(name)]
    if "`suffix'" != "" local suffix = "_`suffix'" 
    if "`if'" != "" local andif = "& (`if')"
    else local andif = ""
    if "`if'" != "" local if = "if (`if')"
    if "`weight'" != "" local wgtexp = "[`weight'=`exp']"
    else local wgtexp = ""
    if "`weight'" != "" local awgt = "[aw=`exp']"
    else local awgt = ""
    tempvar interaction
    tempvar s1_hat
    tempvar s0_hat
    tempvar tau_hat
    tempvar d_help
    tempvar d_hurt
    tempvar d_0
    tempvar tmp_neg
    tempvar pct_utail
    tempvar pct_ltail
    tempvar ll
    tempvar ul
    // Generating survey response variable
    egen `interaction' = group(`tight') `if'
    reg `select' i.(`interaction') `wgtexp' if `treatment' == 1 `andif'
    predict `s1_hat' `if', xb
    reg `select' i.(`interaction') `wgtexp' if `treatment' == 0  `andif'
    predict `s0_hat' `if', xb
    gen `tau_hat' = `s1_hat'-`s0_hat' `if'
    gen `d_help'  = `tau_hat' > 0 `if'
    gen `d_hurt'  = `tau_hat' < 0 `if'
    gen `d_0'     = `tau_hat' == 0 `if'

    sort `tight' `treatment' `varlist'
    di 1
    gen `tmp_neg' = -`varlist' `if'
    bys `tight' `treatment' (`varlist'): cumul `varlist' `awgt' if !mi(`varlist') `andif', equal gen(`pct_ltail')
    bys `tight' `treatment' (`varlist'): cumul `tmp_neg' `awgt' if !mi(`varlist') `andif', equal gen(`pct_utail')
    di 2
    *** Lee upper bound weights
    cap drop wgt_ub`suffix'
    gen wgt_ub`suffix' = 0 `if'
    replace wgt_ub`suffix' = 1 if (`d_0' == 1) `andif'
    ** Upper bound for X(help) - Trimming the lower tail of treatment
    replace wgt_ub`suffix' = 1 if (`select' == 1 & `d_help' == 1 & `treatment' == 1 & `pct_ltail' > abs(`tau_hat')) `andif'
    * Reweighting for non-exact quantiles
    bys `tight' `treatment' (`varlist'): ///
        egen `ll' = ///
            max(0*(`pct_ltail'>=abs(`tau_hat')) + ///
            `pct_ltail'*(`pct_ltail'<abs(`tau_hat'))) ///
        if (`select' == 1 & `d_help' == 1 & `treatment' == 1) `andif'
    bys `tight' `treatment' (`varlist'): ///
        egen `ul' = ///
            min(1*(`pct_ltail'<abs(`tau_hat')) + ///
            `pct_ltail'*(`pct_ltail'>=abs(`tau_hat'))) ///
        if (`select' == 1 & `d_help' == 1 & `treatment' == 1 ) `andif'
    replace wgt_ub`suffix' = (`ul' - abs(`tau_hat')) / (`ul'-`ll') ///
        if (`select' == 1 & `d_help' == 1 & `treatment' == 1 & `pct_ltail' > `ll' & `pct_ltail' <= `ul') `andif'
    drop `ul' `ll'    
    di 3
    ** Upper bound for X(help) - Full Control
    replace wgt_ub`suffix' = 1 if (`select' == 1 & `d_help' == 1 & `treatment' == 0) `andif'

    ** Upper bound for X(hurt) - Full `treatment'
    replace wgt_ub`suffix' = 1 if (`select' == 1 & `d_hurt' == 1 & `treatment' == 1) `andif'
    di 4
    ** Upper bound for X(hurt) - Trimming the upper tail of control
    replace wgt_ub`suffix' = 1 if (`select' == 1 & `d_hurt' == 1 & `treatment' == 0 & `pct_utail' > abs(`tau_hat')) `andif'
    * Reweighting for non-exact quantiles
    bys `tight' `treatment' (`varlist'): ///
        egen `ll' = ///
            max(0*(`pct_utail'>=abs(`tau_hat')) + ///
            `pct_utail'*(`pct_utail'<abs(`tau_hat'))) ///
        if (`select' == 1 & `d_hurt' == 1 & `treatment' == 0) `andif'
    bys `tight' `treatment' (`varlist'): ///
        egen `ul' = ///
            min(1*(`pct_utail'<abs(`tau_hat')) + ///
            `pct_utail'*(`pct_utail'>=abs(`tau_hat'))) ///
        if (`select' == 1 & `d_hurt' == 1 & `treatment' == 0) `andif'
    replace wgt_ub`suffix' = (`ul' - abs(`tau_hat')) / (`ul'-`ll') ///
        if (`select' == 1 & `d_hurt' == 1 & `treatment' == 0 & `pct_utail' > `ll' & `pct_utail' <= `ul') `andif'
    drop `ul' `ll'    

    *** Lee lower bound weights
    cap drop wgt_lb`suffix'
    gen wgt_lb`suffix' = 0 `if'
    replace wgt_lb`suffix' = 1 if (`d_0' == 1) `andif'
    ** Lower bound for X(help) - Trimming the upper tail of treatment
    replace wgt_lb`suffix' = 1 if (`select' == 1 & `d_help' == 1 & `treatment' == 1 & `pct_utail' > abs(`tau_hat')) `andif'
    * Reweighting when there are ties
    bys `tight' `treatment' (`varlist'): ///
        egen `ll' = ///
            max(0*(`pct_utail'>=abs(`tau_hat')) + ///
            `pct_utail'*(`pct_utail'<abs(`tau_hat'))) ///
        if (`select' == 1 & `d_help' == 1 & `treatment' == 1) `andif'
    bys `tight' `treatment' (`varlist'): ///
        egen `ul' = ///
            min(1*(`pct_utail'<abs(`tau_hat')) + ///
            `pct_utail'*(`pct_utail'>=abs(`tau_hat'))) ///
        if (`select' == 1 & `d_help' == 1 & `treatment' == 1) `andif' 
    replace wgt_lb`suffix' = (`ul' - abs(`tau_hat')) / (`ul'-`ll') ///
        if (`select' == 1 & `d_help' == 1 & `treatment' == 1 & `pct_utail' > `ll' & `pct_utail' <= `ul') `andif'
    drop `ul' `ll'     

    ** Lower bound for X(help) - Full Control
    replace wgt_lb`suffix' = 1 if (`select' == 1 & `d_help' == 1 & `treatment' == 0) `andif'

    ** Lower bound for X(hurt) - Full `treatment'
    replace wgt_lb`suffix' = 1 if (`select' == 1 & `d_hurt' == 1 & `treatment' == 1) `andif'

    ** Lower bound for X(hurt) - Trimming the lower tail of control
    replace wgt_lb`suffix' = 1 if (`select' == 1 & `d_hurt' == 1 & `treatment' == 0 & `pct_ltail' > abs(`tau_hat')) `andif'
    * Reweighting for non-exact quantiles
    bys `tight' `treatment' (`varlist'): ///
        egen `ll' = ///
            max(0*(`pct_ltail'>=abs(`tau_hat')) + ///
            `pct_ltail'*(`pct_ltail'<abs(`tau_hat'))) ///
        if (`select' == 1 & `d_hurt' == 1 & `treatment' == 0) `andif'
    bys `tight' `treatment' (`varlist'): ///
        egen `ul' = ///
            min(1*(`pct_ltail'<abs(`tau_hat')) + ///
            `pct_ltail'*(`pct_ltail'>=abs(`tau_hat'))) ///
        if (`select' == 1 & `d_hurt' == 1 & `treatment' == 0) `andif'
    replace wgt_lb`suffix' = (`ul' - abs(`tau_hat')) / (`ul'-`ll') ///
        if (`select' == 1 & `d_hurt' == 1 & `treatment' == 0 & `pct_ltail' > `ll' & `pct_ltail' <= `ul') `andif'
    drop `ul' `ll'  
    drop `tmp_neg' `pct_utail' `pct_ltail' 
    drop `interaction' `s1_hat' `s0_hat' `tau_hat' `d_help' `d_hurt' `d_0'
    replace wgt_ub`suffix' = . if mi(`varlist') `andif'
    replace wgt_lb`suffix' = . if mi(`varlist') `andif'

    di "`weight'"
    if "`weight'" != ""{
        replace wgt_ub`suffix' = wgt_ub`suffix'*`exp'
        replace wgt_lb`suffix' = wgt_lb`suffix'*`exp'
    }
end
