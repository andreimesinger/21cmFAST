pro show_slice, file_name, out_file, redshift, neutral_fraction, dim, color_bar, dt_type, MIN=min_in, MAX=max_in, expand_fact
;
; usage: show_slice, 'filename', 'outfile', 'redshift',
; 'netural_fraction', dim, color_bar, dt_type, [MIN=min_val],
; [MAX=max_val], expand_fact
;
;    file_name, out_file, redshift, neutral_fraction are all strings!!
;
;    color_bar = 0 no color bar
;    color_bar = 1 show color bar
;
;    dt_type = 1 char box
;    dt_type = 4 float box
;
;    expand_fact is the expansion factor, i.e. ratio of final image
;    dimensions to input array dimensions, dim


;  Read in image
box=read_binary(file_name, data_dims=[dim, dim, dim], data_type=dt_type, endian='little')
mid = dim/2 + 0.5
final_size = dim*expand_fact

;box = box*(-1)+1 ;for in_halo
slice_img = extract_slice(box, dim, dim, mid,10,mid, 90,0,0)

min1 = MIN(slice_img)
max1 = MAX(slice_img)

print, 'min, max =', min1, max1

IF KEYWORD_SET(min_in) THEN BEGIN
  min1 = min_in
ENDIF
IF KEYWORD_SET(max_in) THEN BEGIN
  max1 = max_in
ENDIF

data_dim = size(slice_img)
n1 = data_dim(1)
n2 = data_dim(2)

IF (color_bar eq 1) THEN BEGIN
  nadd = ROUND(dim/18)
  n3 = n2+nadd
  data3 = FLTARR(n1, n3)
  data3(*, *) = min1
  print, size(data3)
  data3(0:n1-1, nadd:n3-1) = slice_img

  FOR i=1, n1 DO BEGIN
    data3(i-1, 0:nadd-5) = min1 + (i-0.5)*(max1-min1)/n1
  ENDFOR

  slice_img = data3
  n2 = n3
  
ENDIF

byte_img = BYTSCL(CONGRID (slice_img, final_size, final_size, /INTERP), MIN=min1, MAX=max1, TOP=!D.TABLE_SIZE)

device, decompose=0
;loadct, 0
loadct, 13
WINDOW, 0, retain=2, XSIZE = final_size, YSIZE = final_size
tvscl, byte_img
loadct, 0
redshift_tag='z=' + redshift
nf_tag = 'x!IHI!N=' + neutral_fraction

IF (color_bar eq 1) THEN BEGIN
    xyouts, [final_size/100, final_size/2.1, final_size/1.17], [final_size/100, final_size/100, final_size/100], ['0mK', '50mK', '>100mK'], charsize=0.0033*final_size,color=255, charthick=0.004*final_size, /DEVICE
xyouts, [final_size/1.185, final_size/1.21], [final_size/1.04, final_size/1.09], [redshift_tag, nf_tag], charsize=0.004*final_size,color=255, charthick=0.005*final_size, /DEVICE
ENDIF

void = TVRead(Filename=out_file, /NODIALOG, /JPEG)

END
