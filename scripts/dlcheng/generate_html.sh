SF_DIR=/home/dlcheng/SDSS/G.out
cd $SF_DIR
#rm /home/dlcheng/public_html/images/*
for d in * ; do
    echo "$d"
    if [ -d "$d" ]; 
    then
        convert "$d/$d-D_clusMask.png" -threshold 10% "$d/$d-D3_theshold_clusMask.png"
        #convert "$d/$d-C_preproc.png" "$d/$d-D3_theshold_clusMask.png" "$d/$d-D2_cropMasked.png" +append "$d-CD.png"
        convert "$d/$d-A_input.png" -resize 256x256 "$d/$d-A_input.png"
        convert "$d/$d-L2_residual.png" -resize 256x256 "$d/$d-L2_residual.png"
        convert "$d/$d-L4_clusMask.png" -resize 256x256 "$d/$d-L4_clusMask.png"
        convert "$d/$d-A_input.png" "$d/$d-L2_residual.png" "$d/$d-L4_clusMask.png" +append "/home/dlcheng/public_html/images/$d-CD.png"
        rm "$d/$d-D3_theshold_clusMask.png" -f
    fi
done
chmod 644 ~dlcheng/public_html/images/*
