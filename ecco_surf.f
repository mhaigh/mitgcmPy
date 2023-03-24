           tmpfld=theta(i,j,k,bi,bj)

          tmpmsk=gencost_mskCsurf(i,j,bi,bj,kgen)*
     &               gencost_mskVertical(k,kgen)


            gencost_storefld(i,j,bi,bj,kgen) = 
            gencost_storefld(i,j,bi,bj,kgen)
     &          +tmpmsk*tmpfld*tmpvol
            areavolTile(bi,bj)=areavolTile(bi,bj)
     &          +tmpmsk2*eccoVol_0(i,j,k,bi,bj)



