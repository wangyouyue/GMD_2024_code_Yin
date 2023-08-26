load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRF_contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"

begin
    setfileoption("nc","Format","NetCDF4Classic")

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
    L = 2.5e6
    cp = 1004.0
    Rd = 287.0
    Rv = 461.5
    pref = 1000.0
    g = 9.81

    case_name = "SDM64_50x50x5_mvexchg"
    diri = "../"   ; input directory

    ; Define dimensions.
    dim_names = (/ "x",  "y",  "z", "time" /)
    dim_sizes = (/ I_MAX ,  J_MAX ,  K_MAX ,  T_MAX  /)
    dimUnlim = (/ False , False,  False ,  False  /)

    do fi = 0,MPI - 1
        f = addfile(diri+"history.pe"+ tostring_with_format(fi,"%06d") + ".nc" , "r")
        pres_ave = new((/T_MAX,K_MAX/),"double")
        temp_ave = new((/T_MAX,K_MAX/),"double")
        lwpt = new((/K_MAX,J_MAX,I_MAX/),"double")
        lwpt_h = new((/T_MAX,K_MAX/),"double")
        ql_h = new((/T_MAX,K_MAX/),"double")
        qr_h = new((/T_MAX,K_MAX/),"double")
        qt_h = new((/T_MAX,K_MAX/),"double")
        tke_sgs_h = new((/T_MAX,K_MAX/),"double")
        tke_rs_h = new((/T_MAX,K_MAX/),"double")
        u_ave = new((/T_MAX,K_MAX/),"double")
        v_ave = new((/T_MAX,K_MAX/),"double")
        w_ave = new((/T_MAX,K_MAX/),"double")
        sum_u2 = new((/T_MAX,K_MAX/),"double")
        sum_v2 = new((/T_MAX,K_MAX/),"double")
        sum_w2 = new((/T_MAX,K_MAX/),"double")
        sum_w3 = new((/T_MAX,K_MAX/),"double")
        w_theta_v = new((/T_MAX,K_MAX/),"double")
        dens_ave = new((/T_MAX,K_MAX/),"double")
        pt_ave = new((/T_MAX,K_MAX/),"double")
        theta_v_ave = new((/T_MAX,K_MAX/),"double")
        temp_v_ave = new((/T_MAX,K_MAX/),"double")
        w_h = new((/T_MAX,K_MAX/),"double")
        rh_h = new((/T_MAX,K_MAX/),"double")
        nc_h = new((/T_MAX,K_MAX/),"double")
        cf_h = new((/T_MAX,K_MAX/),"float")
        invz_ave = new((/T_MAX/),"double")
        lwp_t = new((/T_MAX/),"double")
        cf_t = new((/T_MAX/),"float")
        rh_t = new((/T_MAX/),"double")
        ncc = new((/T_MAX/),"integer")
        n_nnt = new((/T_MAX/),"integer")
        nc_t = new((/T_MAX/),"double")
        sfc_precip = new((/T_MAX/),"double")
        n_nnk = new((/T_MAX,K_MAX/),"integer")
        theta_v = new((/T_MAX, K_MAX,J_MAX,I_MAX/),"double")
        ; Define variables.
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
        u_ave@long_name  = "velocity u profile"
        u_ave@units = "m/s"
        v_ave@long_name  = "velocity v profile"
        v_ave@units = "m/s"
        w_ave@long_name  = "velocity w profile"
        w_ave@units = "m/s"
        dens_ave@long_name  = "density profile"
        sum_u2@long_name  = "<u^2> profile"
        sum_u2@units = "m2/s2"
        sum_v2@long_name  = "<v^2> profile"
        sum_v2@units = "m2/s2"
        sum_w2@long_name  = "<w^2> profile"
        sum_w2@units = "m2/s2"
        sum_w3@long_name  = "<w^3> profile"
        sum_w3@units = "m3/s3"
        w_theta_v@long_name  = "<w*thetav> profile"
        w_theta_v@units = "m*k/s"
        dens_ave@units = "kg/m3"
        pt_ave@long_name  = "potential temp. profile"
        pt_ave@units = "K"
        theta_v_ave@long_name  = "virtual potential temp. profile"
        theta_v_ave@units = "K"
        temp_v_ave@long_name  = "virtual temp. profile"
        temp_v_ave@units = "K"
        w_h@long_name  = "w profile"
        w_h@units = "m/s"
        rh_h@long_name  = "relative humidity profile"
        rh_h@units = "%"
        nc_h@long_name  = "cloud droplet number concentration profile in cloudy cells"
        nc_h@units = "cm-3"
        cf_h@long_name  = "cloud fraction profile"
        cf_h@units = ""
        n_nnk@long_name  = "number of grid cells profile"
        n_nnk@units = ""
        invz_ave@long_name  = "inversion height time series"
        invz_ave@units = "m"
        lwp_t@long_name  = "liquid water path time series"
        lwp_t@units = "kg/m2"
        cf_t@long_name  = "cloud cover time series"
        cf_t@units = ""
        rh_t@long_name  = "relative humidity time series"
        rh_t@units = "%"
        ncc@long_name  = "number of cloudy cells series"
        ncc@units = ""
        nc_t@long_name  = "cloud droplet number concentration time series in cloudy cells"
        nc_t@units = "cm-3"
        n_nnt@long_name  = "number of grid cells time series"
        n_nnt@units = ""
        sfc_precip@long_name  = "surface precipitation time series"
        sfc_precip@units = "mm/d"
        sfc_precip = 0
        var_names2D = (/ "pres_ave", "temp_ave", "lwpt_h", "ql_h", "qr_h", "qt_h", "tke_sgs_h", "tke_rs_h", "u_ave", "v_ave", "w_ave", "dens_ave", "pt_ave", "theta_v_ave", "temp_v_ave", "w_h", "rh_h", "sum_u2", "sum_v2", "sum_w2", "sum_w3", "w_theta_v" /)
        varvar_types2D = (/ "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double" /)
        var_names1D = (/ "invz_ave", "lwp_t", "cf_t", "rh_t", "ncc","sfc_precip" /)
        varvar_types1D = (/ "double", "double", "float", "double", "integer", "double" /)
        w_theta_v = 0
        do ti = 0, T_MAX - 1
            PRES = f->PRES(ti,:,:,:)
            TEMP = f->T(ti,:,:,:)  ; temperature, K
            QC_sd   = f->QC_sd(ti,:,:,:)
            QR_sd   = f->QR_sd(ti,:,:,:)
            QLIQ = (QC_sd + QR_sd) * 1000 ; kg/kg -> g/kg
            delete([/QR_sd, QC_sd/])
            QV   = f->QV(ti,:,:,:)  ; kg/kg
            QR   = f->QR_sd(ti,:,:,:) * 1000 ; kg/kg -> g/kg
            QTOT = QLIQ + QV * 1000 ; kg/kg -> g/kg
            TKE_SMG = f->TKE_SMG(ti,:,:,:)   ; turbulent kinetic energy (Smagorinsky), m2/s2
            TKE_RS  = f->TKE_RS(ti,:,:,:)    ; resolved scale TKE, m2/s2
            U       = f->U(ti,:,:,:)
            V       = f->V(ti,:,:,:)
            W       = f->W(ti,:,:,:)
            PT      = f->PT(ti,:,:,:)  ;potential temp.
            DENS    = f->DENS(ti,:,:,:) ;density, kg/m3
            RH      = f->RH(ti,:,:,:)
            RAIN      = f->RAIN(ti,:,:) ; surface rain rate (total), kg/m2/s
            sfc_precip(ti) = sfc_precip(ti) + sum(RAIN) * 86400 ; mm/d
            lwpt = new((/K_MAX,J_MAX,I_MAX/),"double")
            temp_v = new((/K_MAX,J_MAX,I_MAX/),"double")
            do zi = 0, K_MAX - 1
                pres_ave(ti,zi) = avg(PRES(zi,:,:))
                temp_ave(ti,zi) = avg(TEMP(zi,:,:))
                ql_h(ti,zi) = avg(QLIQ(zi,:,:))
                qr_h(ti,zi) = avg(QR(zi,:,:))
                qt_h(ti,zi) = avg(QTOT(zi,:,:))
                tke_sgs_h(ti,zi) = avg(TKE_SMG(zi,:,:))
                tke_rs_h(ti,zi)  = avg(TKE_RS(zi,:,:))
                u_ave(ti,zi) = avg(U(zi,:,:))
                v_ave(ti,zi) = avg(V(zi,:,:))
                w_ave(ti,zi) = avg(W(zi,:,:))
                sum_u2(ti,zi) = sum(U(zi,:,:)^2)
                sum_v2(ti,zi) = sum(V(zi,:,:)^2)
                sum_w2(ti,zi) = sum(W(zi,:,:)^2)
                sum_w3(ti,zi) = sum(W(zi,:,:)^3)
                dens_ave(ti,zi) = avg(DENS(zi,:,:))
                pt_ave(ti,zi) = avg(PT(zi,:,:))
                w_h(ti,zi) = avg(W(zi,:,:))
                rh_h(ti,zi) = avg(RH(zi,:,:))
                do yi = 0, J_MAX - 1
                    do xi = 0, I_MAX - 1
                        lwpt(zi,yi,xi) = ((pref * 100/PRES(zi,yi,xi))^(Rd/cp)) * (TEMP(zi,yi,xi) - L * QLIQ(zi,yi,xi)/(cp * 1000))
                        theta_v(ti,zi,yi,xi) = PT(zi,yi,xi) * (1 + 0.608 * QV(zi,yi,xi) - QLIQ(zi,yi,xi)/1000)
                        temp_v(zi,yi,xi) = TEMP(zi,yi,xi) * (1 + 0.608 * QV(zi,yi,xi))
                        w_theta_v(ti,zi) = w_theta_v(ti,zi) + theta_v(ti,zi,yi,xi) * W(zi,yi,xi)
                    end do
                end do
                lwpt_h(ti,zi) = avg(lwpt(zi,:,:))
                theta_v_ave(ti,zi) = avg(theta_v(ti,zi,:,:))
                temp_v_ave(ti,zi) = avg(temp_v(zi,:,:))
            end do
            delete([/lwpt, temp_v, PRES, TEMP, QLIQ, TKE_SMG, TKE_RS, U, V, W, PT, QV, QR, RAIN/])
            z       = f->z
            QHYD_sd    = f->QHYD_sd(ti,:,:,:) 
            QHYD_sd    = QHYD_sd * 1000
            NC      = f->NC(ti,:,:,:)       ; num/m
            NC = NC*1e-6           ; num/cm3
            invz_ave(ti) = 0
            mass = new((/K_MAX,J_MAX,I_MAX/),"double")
            do yi = 0, J_MAX - 1
                do xi = 0, I_MAX - 1
                    do zi = 0, K_MAX - 1
                        mass(zi,yi,xi) = QHYD_sd(zi,yi,xi)*DENS(zi,yi,xi)*(DX*DY*DZ)*1e-3
                        if (QTOT(zi,yi,xi) .lt. 8.0) then
                            invz_ave(ti) = invz_ave(ti) + z(zi)
                            break
                        end if
                    end do
                end do
            end do
            delete([/QTOT, z, QHYD_sd, DENS/])
            index = ind(ndtooned(NC) .ge. 20)
            if all(ismissing(index)) then
                rh_t(ti) = 0
                ncc(ti) = 0
            else
                RH1D = ndtooned(RH)
                rh_t(ti) = sum(RH1D(index))
                ncc(ti) = dimsizes(index)
                delete(RH1D)
            end if
            delete([/RH, NC, index/])
            cf_t(ti) = num(dim_sum_n(mass, 0)/(DX*DY) .ge. 0.02)
            lwp_t(ti) = sum(mass)/(I_MAX*DX*J_MAX*DY)
            delete(mass)
            invz_ave(ti) = invz_ave(ti)/(I_MAX * J_MAX)
        end do
        delete(theta_v)
        system("rm -f " + case_name + ".pe" + tostring_with_format(fi,"%06d") + ".nc")
        fout = addfile(case_name + ".pe" + tostring_with_format(fi,"%06d") + ".nc","c")
        filedimdef( fout, dim_names, dim_sizes, dimUnlim )
        filevardef(fout, var_names2D, varvar_types2D, (/  "time", "z" /))
        filevardef(fout, var_names1D, varvar_types1D, (/  "time" /))
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
        filevarattdef( fout, "sum_u2", sum_u2 )
        filevarattdef( fout, "sum_v2", sum_v2 )
        filevarattdef( fout, "sum_w2", sum_w2 )
        filevarattdef( fout, "sum_w3", sum_w3 )
        filevarattdef( fout, "w_theta_v", w_theta_v )
        filevarattdef( fout, "dens_ave", dens_ave )
        filevarattdef( fout, "pt_ave", pt_ave )
        filevarattdef( fout, "theta_v_ave", theta_v_ave )
        filevarattdef( fout, "temp_v_ave", temp_v_ave )
        filevarattdef( fout, "w_h", w_h )
        filevarattdef( fout, "rh_h", rh_h )
        filevarattdef( fout, "invz_ave", invz_ave )
        filevarattdef( fout, "lwp_t", lwp_t )
        filevarattdef( fout, "cf_t", cf_t )
        filevarattdef( fout, "rh_t", rh_t )
        filevarattdef( fout, "sfc_precip", sfc_precip )
        filevarattdef( fout, "ncc", ncc )

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
        fout->sum_u2    = (/sum_u2/)
        fout->sum_v2    = (/sum_v2/)
        fout->sum_w2    = (/sum_w2/)
        fout->sum_w3    = (/sum_w3/)
        fout->w_theta_v = (/w_theta_v/)
        fout->dens_ave  = (/dens_ave/)
        fout->pt_ave    = (/pt_ave/)
        fout->theta_v_ave    = (/theta_v_ave/)
        fout->temp_v_ave     = (/temp_v_ave/)
        fout->w_h       = (/w_h/)
        fout->rh_h      = (/rh_h/)
        fout->invz_ave  = (/invz_ave/)
        fout->lwp_t     = (/lwp_t/)
        fout->cf_t      = (/cf_t/)
        fout->rh_t      = (/rh_t/)
        fout->sfc_precip      = (/sfc_precip/)
        fout->ncc       = (/ncc/)
        delete([/ fout, sfc_precip, pres_ave, temp_ave, lwpt_h, ql_h, qr_h, qt_h, tke_sgs_h, tke_rs_h, invz_ave, u_ave, v_ave, w_ave, dens_ave, pt_ave, theta_v_ave, temp_v_ave, lwp_t, cf_t, rh_h, w_h, rh_t, ncc, w_theta_v, sum_u2, sum_v2, sum_w2, sum_w3/])
        ; 3D (t, y, x)
        QHYD_sd = f->QHYD_sd
        QHYD_sd = QHYD_sd * 1000
        DENS = f->DENS    ;density, kg/m3
        mass_2d = new((/T_MAX, K_MAX, J_MAX, I_MAX/),"double")
        mass_2d = QHYD_sd*DENS*(DX*DY*DZ)*1e-3
        delete([/QHYD_sd, DENS/])
        lwp_2d = new((/T_MAX, J_MAX, I_MAX/),"double")
        lwp_2d = dim_sum_n(mass_2d, 1)/(DX*DY)
        lwp_2d@long_name  = "2D liquid water path"
        lwp_2d@units = "kg/m2"
        var_names3D = (/ "lwp_2d" /)
        varvar_types3D = (/ "double" /)
        fout = addfile(case_name + ".pe" + tostring_with_format(fi,"%06d") + ".nc","w")
        filevardef(fout, var_names3D, varvar_types3D, (/  "time", "y", "x" /))
        filevarattdef( fout, "lwp_2d", lwp_2d )
        fout->lwp_2d      = (/lwp_2d/)
        delete([/fout, lwp_2d/])
        ; 4D (t, z, y, x)
        NC   = f->NC        ; num/kg
        delete(f)
        NC = NC * 1e-6  ; num/cm3
        cloud = new((/T_MAX, K_MAX, J_MAX, I_MAX/),"integer")
        cloud = 0
        NC1D = ndtooned(NC)
        index = ind(NC1D .ge. 20)
        delete(NC1D)
        cf_h = 0
        nc_h = 0
        n_nnk = 0
        nc_t = 0
        n_nnt = 0
        if (all(.not.ismissing(index))) then
            indices  = ind_resolve(index, dimsizes(NC))
            dim_ida  = dimsizes(indices)
            npts     = dim_ida(0)       ; number of elements >= 20
            ndim     = dim_ida(1)       ; rank of NC
            do n=0, npts-1
                cloud(indices(n,0), indices(n,1), indices(n,2), indices(n,3)) = 1
            end do
            start = 0
            do ti = 0, T_MAX - 1
                nnt = ind((indices(start:,0) .eq. ti))
                if (all(.not.ismissing(nnt))) then
                    n_nnt(ti) = dimsizes(nnt)
                    ll = nnt(n_nnt(ti)-1)
                    nnt = nnt + start
                    start = start +  ll + 1
                    if (start .ge. npts) then
                        start = npts - 1
                    end if
                    do ni = 0, n_nnt(ti) - 1
                        nc_t(ti) = nc_t(ti) + NC(ti, indices(nnt(ni),1), indices(nnt(ni),2), indices(nnt(ni),3))
                    end do
                    do zi = 0, K_MAX - 1
                        cf_h(ti, zi) = sum(cloud(ti, zi, :, :))
                        nnk = ind((indices(nnt,1) .eq. zi))
                        if (all(.not.ismissing(nnk))) then
                            n_nnk(ti, zi) = dimsizes(nnk)
                            do nni = 0, n_nnk(ti, zi) - 1
                                nc_h(ti, zi) = nc_h(ti, zi) + NC(ti, zi, indices(nnt(nnk(nni)),2), indices(nnt(nnk(nni)),3))
                            end do
                        end if
                        delete(nnk)
                    end do
                else
                    do zi = 0, K_MAX - 1
                        cf_h(ti, zi) = sum(cloud(ti, zi, :, :))
                    end do
                end if
                delete(nnt)
            end do
            delete([/indices, dim_ida, npts, ndim/])
        end if
        delete(index)
        fout = addfile(case_name + ".pe" + tostring_with_format(fi,"%06d") + ".nc","w")
        filevardef(fout, (/"nc_t", "n_nnt"/), (/ "double", "integer"/), (/  "time" /))
        filevarattdef( fout, "nc_t", nc_t )
        filevarattdef( fout, "n_nnt", n_nnt )
        fout->nc_t        = (/nc_t/)
        fout->n_nnt       = (/n_nnt/)
        delete([/nc_t, n_nnt/])
        filevardef(fout, (/"cf_h", "nc_h", "n_nnk"/), (/ "float", "double", "integer"/), (/  "time", "z" /))
        filevarattdef( fout, "cf_h", cf_h )
        filevarattdef( fout, "nc_h", nc_h )
        filevarattdef( fout, "n_nnk", n_nnk )
        fout->cf_h       = (/cf_h/)
        fout->nc_h       = (/nc_h/)
        fout->n_nnk       = (/n_nnk/)
        delete([/nc_h, n_nnk, cf_h/])
        NC@long_name  = "3D cloud droplet number concentration"
        NC@units = "num/cm3"
        cloud@long_name  = "cloudy cells: 1, non-cloudy cells: 0"
        cloud@units = ""
        var_names4D = (/ "NC", "cloud" /)
        varvar_types4D = (/ "double", "integer" /)
        filevardef(fout, var_names4D, varvar_types4D, (/  "time", "z", "y", "x" /))
        filevarattdef( fout, "NC", NC )
        filevarattdef( fout, "cloud", cloud )
        fout->NC        = (/NC/)
        fout->cloud     = (/cloud/)
        delete(NC)
        ; cloud top/base height
        cloud_base = new((/T_MAX,J_MAX,I_MAX/),"float")
        cloud_base = 0
        cloud_base@long_name  = "cloud base height time series"
        cloud_base@units = "m"
        cloud_top = new((/T_MAX,J_MAX,I_MAX/),"float")
        cloud_top = 0
        cloud_top@long_name  = "cloud top height time series"
        cloud_top@units = "m"
        n_col = new((/T_MAX/),"integer")
        n_col = 0
        n_col@long_name  = "number of cloudy columns time series"
        n_col@units = ""
        do ti = 0, T_MAX - 1
            do yi = 0, J_MAX - 1
                do xi = 0, I_MAX - 1
                    if (sum(cloud(ti,:,yi,xi)) .gt. 0) then 
                        n_col(ti) = n_col(ti) + 1
                        zi = 0
                        do while(cloud(ti,zi,yi,xi) .eq. 0)
                            zi = zi + 1
                        end do
                        cloud_base(ti,yi,xi) = zi * DZ - DZ/2.0
                        zi = K_MAX - 1
                        do while(cloud(ti,zi,yi,xi) .eq. 0)
                            zi = zi - 1
                        end do
                        cloud_top(ti,yi,xi) = zi * DZ - DZ/2.0
                    end if
                end do
            end do
        end do
        delete(cloud)
        filevardef(fout, (/"n_col"/), (/ "integer"/), (/  "time" /))
        filevarattdef( fout, "n_col", n_col )
        fout->n_col        = (/n_col/)
        delete(n_col)
        filevardef(fout, (/"cloud_base", "cloud_top"/), (/ "float", "float"/), (/  "time", "y", "x" /))
        filevarattdef( fout, "cloud_base", cloud_base )
        filevarattdef( fout, "cloud_top", cloud_top )
        fout->cloud_top        = (/cloud_top/)
        fout->cloud_base       = (/cloud_base/)
        delete([/fout, cloud_base, cloud_top/])
    end do
end
