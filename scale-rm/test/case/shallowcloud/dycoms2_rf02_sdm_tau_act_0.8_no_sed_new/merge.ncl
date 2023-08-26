load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRF_contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"

begin
	setfileoption("nc","Format","NetCDF4Classic")
    ntime = 360
    PRC_NUM_X = 20
    PRC_NUM_Y = 20
    MPI = PRC_NUM_X*PRC_NUM_Y
    I_MAX = 6
    J_MAX = 6
    K_MAX = 300
    T_MAX = 361
    DX = 50
    DY = 50
    DZ = 5
    DT = 60
    D = 3.75e-6
    g = 9.81

	diri = "./"   ; input directory
    case_name = "SDM64_50x50x5_mvexchg"

    ; Define dimensions.
    dim_names = (/ "x",  "y",  "z", "time", "X", "Y", "ntime"/)
    dim_sizes = (/ I_MAX, J_MAX, K_MAX, T_MAX, I_MAX*PRC_NUM_X, J_MAX*PRC_NUM_Y, ntime+1 /)
    dimUnlim = (/ False, False, False, False, False, False, False /)
    ; 2D
    pres_ave = new((/T_MAX,K_MAX/),"double")
    pres_ave = 0
    temp_ave = new((/T_MAX,K_MAX/),"double")
    temp_ave = 0
    lwpt_h = new((/T_MAX,K_MAX/),"double")
    lwpt_h = 0
    ql_h = new((/T_MAX,K_MAX/),"double")
    ql_h = 0
    qr_h = new((/T_MAX,K_MAX/),"double")
    qr_h = 0
    qt_h = new((/T_MAX,K_MAX/),"double")
    qt_h = 0
    tke_sgs_h = new((/T_MAX,K_MAX/),"double")
    tke_sgs_h = 0
    tke_rs_h = new((/T_MAX,K_MAX/),"double")
    tke_rs_h = 0
	dens_ave = new((/T_MAX,K_MAX/),"double")
	dens_ave = 0
    pt_ave = new((/T_MAX,K_MAX/),"double")
    pt_ave = 0
    rh_h = new((/T_MAX,K_MAX/),"double")
    rh_h = 0
    theta_v_ave = new((/T_MAX,K_MAX/),"double")
    theta_v_ave = 0
    temp_v_ave = new((/T_MAX,K_MAX/),"double")
    temp_v_ave = 0
    u_ave = new((/T_MAX,K_MAX/),"double")
    u_ave = 0
    v_ave = new((/T_MAX,K_MAX/),"double")
    v_ave = 0
    w_ave = new((/T_MAX,K_MAX/),"double")
    w_ave = 0
    nc_h = new((/T_MAX,K_MAX/),"double")
    nc_h = 0
    cf_h = new((/T_MAX,K_MAX/),"float")
    cf_h = 0
    u2_ave = new((/T_MAX,K_MAX/),"double")
    u2_ave = 0
    v2_ave = new((/T_MAX,K_MAX/),"double")
    v2_ave = 0
    w2_ave = new((/T_MAX,K_MAX/),"double")
    w2_ave = 0
    w3_ave = new((/T_MAX,K_MAX/),"double")
    w3_ave = 0
    pres_ave@long_name  = "pressure profile"
    pres_ave@units = "Pa"
    temp_ave@long_name  = "temperature profile"
    temp_ave@units = "K"
    lwpt_h@long_name  = "liquid water potential temperature profile"
    lwpt_h@units = "K"
    ql_h@long_name  = "liquid water mixing ratio profile"
    ql_h@units = "g/kg"
    qr_h@long_name  = "rain water mixing ratio profile"
    qr_h@units = "g/kg"
    qt_h@long_name  = "total water mixing ratio profile"
    qt_h@units = "g/kg"
    tke_sgs_h@long_name  = "turbulent kinetic energy (Smagorinsky) profile"
    tke_sgs_h@units = "m2/s2"
    tke_rs_h@long_name  = "resolved scale TKE profile"
    tke_rs_h@units = "m2/s2"
    dens_ave@long_name  = "density profile"
    dens_ave@units = "kg/m3"
    pt_ave@long_name  = "potential temp. profile"
    pt_ave@units = "K"
    rh_h@long_name  = "relative humidity profile"
    rh_h@units = "%"
    theta_v_ave@long_name  = "virtual potential temp. profile"
    theta_v_ave@units = "K"
    temp_v_ave@long_name  = "virtual temp. profile"
    temp_v_ave@units = "K"
    u_ave@long_name  = "velocity u profile"
    u_ave@units = "m/s"
    v_ave@long_name  = "velocity v profile"
    v_ave@units = "m/s"
    w_ave@long_name  = "velocity w profile"
    w_ave@units = "m/s"
    nc_h@long_name  = "cloud droplet number concentration profile in cloudy cells"
    nc_h@units = "cm-3"
    cf_h@long_name  = "cloud fraction profile"
    cf_h@units = ""
    u_val = new((/T_MAX,K_MAX/),"double")
    u_val = 0
    v_val = new((/T_MAX,K_MAX/),"double")
    v_val = 0
    w_val = new((/T_MAX,K_MAX/),"double")
    w_val = 0
    w_sk = new((/T_MAX,K_MAX/),"double")
    w_sk = 0
    u_val@long_name  = "u variance profile"
    u_val@units = "m2/s2"
    v_val@long_name  = "v variance profile"
    v_val@units = "m2/s2"
    w_val@long_name  = "w variance profile"
    w_val@units = "m2/s2"
    w_sk@long_name  = "w skewness profile"
    w_sk@units = "m3/s3"
    tot_boy = new((/T_MAX,K_MAX/),"double")
    tot_boy = 0
    tot_boy@long_name  = "buoyancy production profile"
    tot_boy@units = "cm2/s3"
    var_names2D = (/ "pres_ave", "temp_ave", "lwpt_h", "ql_h", "qr_h", "qt_h", "tke_sgs_h", "tke_rs_h", "u_ave", "v_ave", "w_ave", "dens_ave", "pt_ave", "theta_v_ave", "temp_v_ave", "rh_h", "nc_h", "cf_h", "u_val", "v_val", "w_val", "w_sk", "tot_boy" /)
    varvar_types2D = (/ "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "float", "double", "double", "double", "double", "double" /)
    ; 1D
    invz_ave = new((/T_MAX/),"double")
    invz_ave = 0
    invz_ave@long_name  = "inversion height time series"
    invz_ave@units = "m"
    er = new((/T_MAX/),"double")
    er = 0
    er@long_name  = "entrainment rate time series"
    er@units = "cm/s"
    lwp_t = new((/T_MAX/),"double")
    lwp_t = 0
    lwp_t@long_name  = "liquid water path time series"
    lwp_t@units = "kg/m2"
    rh_t = new((/T_MAX/),"double")
    rh_t = 0
    rh_t@long_name  = "relative humidity time series"
    rh_t@units = "%"
    sfc_precip = new((/T_MAX/),"double")
    sfc_precip = 0
    sfc_precip@long_name  = "surface precipitation time series"
    sfc_precip@units = "mm/d"
    cf_t = new((/T_MAX/),"float")
    cf_t = 0
    cf_t@long_name  = "cloud cover time series"
    cf_t@units = ""
    cloud_top = new((/T_MAX/),"float")
    cloud_top = 0
    cloud_top@long_name  = "cloud top height time series"
    cloud_top@units = "m"
    cloud_base = new((/T_MAX/),"float")
    cloud_base = 0
    cloud_base@long_name  = "cloud base height time series"
    cloud_base@units = "m"
    nc_t = new((/T_MAX/),"double")
    nc_t = 0
    nc_t@long_name  = "cloud droplet number concentration time series in cloudy cells"
    nc_t@units = "cm-3"
    tke_t = new((/T_MAX/),"double")
    tke_t = 0
    tke_t@long_name  = "turbulent kinetic energy time series"
    tke_t@units = "m2/s2"
    w_max = new((/T_MAX/),"double")
    w_max@long_name  = "max. of w variance time series"
    w_max@units = "m2/s2"
    var_names1D = (/ "invz_ave", "er", "lwp_t", "rh_t", "cf_t", "nc_t", "tke_t", "w_max", "sfc_precip", "cloud_top", "cloud_base" /)
    varvar_types1D = (/ "double", "double", "double", "double", "float", "double", "double", "double", "double", "float", "float" /)
    n_nnk = new((/T_MAX,K_MAX/),"integer")
    n_nnk = 0
    ncc = new((/T_MAX/),"integer")
    ncc = 0
    n_nnt = new((/T_MAX/),"integer")
    n_nnt = 0
    n_col = new((/T_MAX/),"integer")
    n_col = 0
    w_theta_v_ave = new((/T_MAX,K_MAX/),"double")
    w_theta_v_ave = 0
    do fi = 0, MPI - 1
    	f = addfile(case_name + ".pe" + tostring_with_format(fi,"%06d") + ".nc", "r")
    	; 2D
    	PRES = f->pres_ave
    	TEMP = f->temp_ave
    	LWPT = f->lwpt_h
    	QL   = f->ql_h
        QR   = f->qr_h
    	QT   = f->qt_h
    	TKE_SGS   = f->tke_sgs_h
    	TKE_RS    = f->tke_rs_h
    	DENS = f->dens_ave
    	PT   = f->pt_ave
    	RH   = f->rh_h
    	theta_v = f->theta_v_ave
    	temp_v  = f->temp_v_ave
    	sum_u    = f->u_ave
    	sum_v    = f->v_ave
    	sum_w    = f->w_ave
        sum_u2    = f->sum_u2
        sum_v2    = f->sum_v2
        sum_w2    = f->sum_w2
        sum_w3    = f->sum_w3
        w_theta_v = f->w_theta_v
    	nc_profile = f->nc_h
    	NNK  = f->n_nnk
    	cf_profile = f->cf_h

    	pres_ave = pres_ave + PRES
    	temp_ave = temp_ave + TEMP
    	lwpt_h = lwpt_h + LWPT
    	ql_h = ql_h + QL
        qr_h = qr_h + QR
    	qt_h = qt_h + QT
    	tke_sgs_h = tke_sgs_h + TKE_SGS
    	tke_rs_h  = tke_rs_h + TKE_RS
    	dens_ave  = dens_ave + DENS
    	pt_ave    = pt_ave + PT
    	rh_h = rh_h + RH
    	theta_v_ave = theta_v_ave + theta_v
    	temp_v_ave  = temp_v_ave + temp_v
    	u_ave = u_ave + sum_u
    	v_ave = v_ave + sum_v
    	w_ave = w_ave + sum_w
        u2_ave = u2_ave + sum_u2
        v2_ave = v2_ave + sum_v2
        w2_ave = w2_ave + sum_w2
        w3_ave = w3_ave + sum_w3
        w_theta_v_ave = w_theta_v_ave + w_theta_v
    	nc_h  = nc_h + nc_profile
    	n_nnk = n_nnk + NNK
    	cf_h  = cf_h + cf_profile
    	delete([/ PRES, TEMP, LWPT, QL, QT, TKE_SGS, TKE_RS, DENS, PT, RH, theta_v, temp_v, sum_u, sum_v, sum_w, nc_profile, NNK, cf_profile, sum_u2, sum_v2, sum_w2, sum_w3, w_theta_v/])
        
        ; 1D
        INVZ = f->invz_ave
        LWP  = f->lwp_t
        RH_TIME= f->rh_t
        RAIN = f->sfc_precip
        NUMBER = f->ncc
        CF   = f->cf_t
        NC   = f->nc_t
        NNT  = f->n_nnt
        col= f->n_col
        top  = f->cloud_top
        base = f->cloud_base
        delete(f)
        invz_ave = invz_ave + INVZ
        lwp_t    = lwp_t + LWP
        rh_t     = rh_t + RH_TIME
        ncc      = ncc + NUMBER
        cf_t     = cf_t + CF
        nc_t     = nc_t + NC
        n_nnt    = n_nnt + NNT
        sfc_precip = sfc_precip + RAIN
        n_col    = n_col + col
        do ti = 1, T_MAX - 1
            cloud_top(ti) = cloud_top(ti) + sum(top(ti,:,:))
            cloud_base(ti) = cloud_base(ti) + sum(base(ti,:,:))
        end do

        delete([/INVZ, LWP, RH_TIME, NUMBER, CF, NC, NNT, RAIN, top, base/])
    end do
    ; 2D
    pres_ave = pres_ave/MPI
    temp_ave = temp_ave/MPI
    lwpt_h   = lwpt_h/MPI
    ql_h     = ql_h/MPI
    qr_h     = qr_h/MPI
    qt_h     = qt_h/MPI
    tke_sgs_h= tke_sgs_h/MPI
    tke_rs_h = tke_rs_h/MPI
    dens_ave = dens_ave/MPI
    pt_ave   = pt_ave/MPI
    rh_h     = rh_h/MPI
    theta_v_ave = theta_v_ave/MPI
    temp_v_ave  = temp_v_ave/MPI
    u_ave    = u_ave/MPI
    v_ave    = v_ave/MPI
    w_ave    = w_ave/MPI
    u2_ave = u2_ave/(MPI*I_MAX*J_MAX)
    v2_ave = v2_ave/(MPI*I_MAX*J_MAX)
    w2_ave = w2_ave/(MPI*I_MAX*J_MAX)
    w3_ave = w3_ave/(MPI*I_MAX*J_MAX)
    w_theta_v_ave = w_theta_v_ave/(MPI*I_MAX*J_MAX)
    cf_h     = cf_h/(MPI*I_MAX*J_MAX)
    u_val    = u2_ave - u_ave^2
    v_val    = v2_ave - v_ave^2
    w_val    = w2_ave - w_ave^2
    w_sk     = w3_ave - w_ave^3 + 3 * w_ave * (w_ave^2 - w2_ave)
    tot_boy = (w_theta_v_ave - w_ave * theta_v_ave)*g*10000.0/temp_v_ave
    delete([/u2_ave, v2_ave, w2_ave, w3_ave, w_theta_v_ave/])
    do ti = 1, T_MAX - 1
    	do zi = 0, K_MAX - 1
    		if (n_nnk(ti, zi) .eq. 0) then
    			nc_h(ti, zi) = 0
    		else
    			nc_h(ti, zi) = nc_h(ti, zi)/n_nnk(ti, zi)
    		end if
    	end do
    end do
    system("rm -f " + case_name + ".nc")
    fout = addfile(case_name + ".nc","c")
    filedimdef( fout, dim_names, dim_sizes, dimUnlim )
    filevardef(fout, var_names2D, varvar_types2D, (/  "time", "z" /))
    filevarattdef( fout, "pres_ave", pres_ave )
    filevarattdef( fout, "temp_ave", temp_ave )
    filevarattdef( fout, "lwpt_h", lwpt_h )
    filevarattdef( fout, "ql_h", ql_h )
    filevarattdef( fout, "qr_h", qr_h )
    filevarattdef( fout, "qt_h", qt_h )
    filevarattdef( fout, "tke_sgs_h", tke_sgs_h )
    filevarattdef( fout, "tke_rs_h", tke_rs_h )
    filevarattdef( fout, "u_ave", u_ave )
    filevarattdef( fout, "v_ave", v_ave )
    filevarattdef( fout, "w_ave", w_ave )
    filevarattdef( fout, "u_val", u_val )
    filevarattdef( fout, "v_val", v_val )
    filevarattdef( fout, "w_val", w_val )
    filevarattdef( fout, "w_sk", w_sk )
    filevarattdef( fout, "dens_ave", dens_ave )
    filevarattdef( fout, "pt_ave", pt_ave )
    filevarattdef( fout, "theta_v_ave", theta_v_ave )
    filevarattdef( fout, "temp_v_ave", temp_v_ave )
    filevarattdef( fout, "rh_h", rh_h )
    filevarattdef( fout, "nc_h", nc_h )
    filevarattdef( fout, "cf_h", cf_h )
    filevarattdef( fout, "tot_boy", tot_boy )

    fout->pres_ave  = (/pres_ave/)
    fout->temp_ave  = (/temp_ave/)
    fout->lwpt_h    = (/lwpt_h/)
    fout->ql_h      = (/ql_h/)
    fout->qr_h      = (/qr_h/)
    fout->qt_h      = (/qt_h/)
    fout->tke_sgs_h = (/tke_sgs_h/)
    fout->tke_rs_h  = (/tke_rs_h/)
    fout->u_ave     = (/u_ave/)
    fout->v_ave     = (/v_ave/)
    fout->w_ave     = (/w_ave/)
    fout->u_val     = (/u_val/)
    fout->v_val     = (/v_val/)
    fout->w_val     = (/w_val/)
    fout->w_sk      = (/w_sk/)
    fout->dens_ave  = (/dens_ave/)
    fout->pt_ave    = (/pt_ave/)
    fout->theta_v_ave    = (/theta_v_ave/)
    fout->temp_v_ave     = (/temp_v_ave/)
    fout->rh_h      = (/rh_h/)
    fout->nc_h      = (/nc_h/)
    fout->cf_h      = (/cf_h/)
    fout->tot_boy   = (/tot_boy/)
    delete([/pres_ave, temp_ave, lwpt_h, ql_h, qr_h, qt_h, tke_sgs_h, tke_rs_h, dens_ave, pt_ave, rh_h, nc_h, cf_h, n_nnk, w_sk, temp_v_ave, tot_boy/])
    ; 1D
    lwp_t = lwp_t/MPI
    invz_ave = invz_ave/MPI
    cf_t = cf_t/(MPI*I_MAX*J_MAX)
    sfc_precip = sfc_precip/(MPI*I_MAX*J_MAX)
    do ti = 1, T_MAX - 1
        if (ncc(ti) .eq. 0) then
            rh_t(ti) = rh_t@_FillValue
        else
            rh_t(ti) = rh_t(ti)/ncc(ti)
        end if
        if (n_nnt(ti) .eq. 0) then
            nc_t(ti) = nc_t@_FillValue
        else
            nc_t(ti) = nc_t(ti)/n_nnt(ti)
        end if
        if (n_col(ti) .eq. 0) then
            cloud_top(ti) = cloud_top@_FillValue
            cloud_base(ti) = cloud_base@_FillValue
        else
            cloud_top(ti) = cloud_top(ti)/n_col(ti)
            cloud_base(ti) = cloud_base(ti)/n_col(ti)
        end if
    end do
    delete(ncc)
    delete(n_col)
    do ti = 1, T_MAX - 1
       er(ti) = (invz_ave(ti) - invz_ave(ti-1))/DT + D * invz_ave(ti)
    end do
	do ti = 1, T_MAX - 1
       w_max(ti) = max(w_val(ti,:))
       do zi = 0, K_MAX - 1
          tke_t(ti) = tke_t(ti) + ((u_val(ti,zi) + v_val(ti,zi) + w_val(ti,zi))/2)*DZ
       end do
	end do
	er = er * 100
	filevardef(fout, var_names1D, varvar_types1D, (/  "time" /))
	filevarattdef( fout, "invz_ave", invz_ave )
	filevarattdef( fout, "er", er )
	filevarattdef( fout, "lwp_t", lwp_t )
	filevarattdef( fout, "rh_t", rh_t )
    filevarattdef( fout, "sfc_precip", sfc_precip )
	filevarattdef( fout, "cf_t", cf_t )
	filevarattdef( fout, "nc_t", nc_t )
    filevarattdef( fout, "tke_t", tke_t )
    filevarattdef( fout, "w_max", w_max )
    filevarattdef( fout, "cloud_top", cloud_top )
    filevarattdef( fout, "cloud_base", cloud_base )

	fout->invz_ave= (/invz_ave/)
	fout->er      = (/er/)
	fout->lwp_t   = (/lwp_t/)
	fout->rh_t    = (/rh_t/)
    fout->sfc_precip    = (/sfc_precip/)
	fout->cf_t    = (/cf_t/)
	fout->nc_t    = (/nc_t/)
    fout->tke_t   = (/tke_t/)
    fout->w_max   = (/w_max/)
    fout->cloud_top   = (/cloud_top/)
    fout->cloud_base   = (/cloud_base/)
	delete([/invz_ave, er, lwp_t, rh_t, cf_t, nc_t, n_nnt, u_val, v_val, w_val, w_max, sfc_precip, tke_t, cloud_base, cloud_top/])
    ; 4D
    files = systemfunc ("ls " + case_name + ".pe*") 
    f = addfiles (files, "r")   
    ListSetType (f, "join")
    lwp_2d = new((/ntime + 1, J_MAX*PRC_NUM_Y,I_MAX*PRC_NUM_X/),"double")
    lwp_2d@long_name  = "liquid water path horizontal distribution"
    lwp_2d@units = "kg/m2"
    do ti = 1, T_MAX - 1, (T_MAX - 1)/ntime
        lwp = f[:]->lwp_2d(:,ti,:,:)
        do fi = 0, MPI - 1
            do xi = 0, I_MAX - 1
                xp = mod(fi,PRC_NUM_X)*I_MAX + xi
                do yi = 0, J_MAX - 1
                    yp = (fi/PRC_NUM_X)*J_MAX + yi
                    lwp_2d(ti,yp,xp) = lwp(fi,yi,xi) * 1
                end do
            end do
        end do
        delete(lwp)
    end do
    delete(f)
    filevardef(fout, "lwp_2d", "double", (/ "ntime", "Y", "X" /))
    filevarattdef( fout, "lwp_2d", lwp_2d )
    fout->lwp_2d  = (/lwp_2d/)
    delete([/fout, lwp_2d/])
end
