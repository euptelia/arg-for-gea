sleep 3h
echo 6x continuous selection test starts
date
SECONDS=0
slim /home/anadem/github/arg-for-gea/slim/second_revision/test/contiuous_nonWF_M3a_glacialHistory_clineMap_recurrentChange_adjustMut.slim
echo 6x continuous selection test ends
date
echo $(($SECONDS/60)) minutes
echo $(($SECONDS/3600)) hours

echo discrete selection test starts
date
SECONDS=0
slim /home/anadem/github/arg-for-gea/slim/second_revision/test/simpleModelForTesting_largeN.slim
echo discrete selection test ends
date
echo $(($SECONDS/60)) minutes
echo $(($SECONDS/3600)) hours


echo discrete neutral test starts
date
SECONDS=0
slim /home/anadem/github/arg-for-gea/slim/second_revision/test/simpleModelForTesting_largeN_neutral.slim
echo discrete neutral test ends
date
echo $(($SECONDS/60)) minutes
echo $(($SECONDS/3600)) hours

echo continuous neutral test starts
date
SECONDS=0
slim /home/anadem/github/arg-for-gea/slim/second_revision/test/contiuous_nonWF_M3a_glacialHistory_clineMap_recurrentChange_adjustMut_neutral.slim
echo continuous neutral test ends
date
echo $(($SECONDS/60)) minutes
echo $(($SECONDS/3600)) hours




