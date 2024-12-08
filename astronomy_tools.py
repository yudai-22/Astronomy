def v2ch(v, w): # v(km/s)をchに変える
    x_tempo, y_tempo, v_tempo   = w.wcs_pix2world(0, 0, 0, 0)
    x_ch, y_ch, v_ch   = w.wcs_world2pix(x_tempo, y_tempo, v*1000.0, 0)
    v_ch = int(round(float(v_ch), 0))
    return v_ch


def ch2v(ch, w):  # チャンネル番号を視線速度 (v: km/s) に変換
    x_tempo, y_tempo, v_tempo = w.wcs_pix2world(0, 0, ch, 0)
    v_kms = v_tempo / 1000.0  # m/s を km/s に変換
    return v_kms


def icrs_to_galactic(ra, dec):#赤経赤緯を銀経銀緯へ
    coord_icrs = SkyCoord(ra=ra * u.deg, dec=dec * u.deg, frame='icrs')  # ICRS座標系(J2000)
    galactic = coord_icrs.galactic
    return galactic.l.deg, galactic.b.deg


def del_header_key(header, keys): # headerのkeyを消す
    import copy
    h = copy.deepcopy(header)
    for k in keys:
        try:
            del h[k]
        except:
            pass
    return h


def make_new_hdu_integ(hdu, v_start_ch, v_end_ch, w): # 指定速度積分強度のhduを作る
    data = hdu.data
    header = hdu.header
    new_data = np.nansum(data[v_start_ch:v_end_ch+1], axis=0)*np.abs(header["CDELT3"])/1000.0
    header = del_header_key(header, ["CRVAL3", "CRPIX3", "CRVAL3", "CDELT3", "CUNIT3", "CTYPE3", "CROTA3", "NAXIS3", "PC1_3", "PC2_3", "PC3_3", "PC3_1", "PC3_2"])
    header["NAXIS"] = 2
    new_hdu = fits.PrimaryHDU(new_data, header)
    return new_hdu


def plot_selected_channel(data, start_ch=None, end_ch=None):#data_shape=(depth, width, height), plot mean intensity
    plt.figure(figsize=(6, 6))
    mean_data = np.nanmean(data, axis=(1, 2))
    
    plt.plot(np.arange(len(mean_data)), mean_data, "k")
    plt.xlabel("Channel")
    plt.xticks(np.arange(0, len(mean_data), 50), fontsize=10)
    plt.ylabel("Mean Intensity [K]")

    if start_ch is not None:
        plt.axvline(x=start_ch, color='red', linestyle='--', label=f'channel = {start_ch}')
    if end_ch is not None:
        plt.axvline(x=end_ch, color='b', linestyle='--', label=f'channel = {end_ch}')

    plt.grid()
    plt.legend()
    plt.title("", fontsize=14)

    plt.show()

    
    
    
