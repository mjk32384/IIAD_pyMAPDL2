## 파일 설명
``functions.py``: tensile test 및 geomtry file 생성.  
``test/test_*.py``: ``functions.py``에 적절한 파라미터를 넘겨주어서, 원하는 시뮬레이션.  
``result_visulization.ipynb``: 시뮬레이션 결과 시각화. 다만 mapdl 인스턴스를 불러온 결과는 다른 컴퓨터에서는 안 보임.

``result/``: 시뮬레이션 결과가 저장되는 폴더. 용량 문제로 github에 업로드는 안하고 있음.  
``Examples/``: pyMAPDL 예제 있는 폴더. 무시해도 됨.  
``temp/``: 이것저것 시도해본 파일. 무시해도 됨.

## 시뮬레이션 방법
Multilinear Isotropic hardening model (MISO Model)을 통해 non-linear 부분을 시뮬레이션 하고, Maximum Normal Stress Theory를 이용해서 fracture 여부를 판단하고 있음.

[신도리코](https://sindoh4u.com/goods/view?no=46&srsltid=AfmBOopVEhJF0xLJstQIMS8SeHpXZjQAaL1keqCtixq9tZ3yZEf-YIp9)에서는 PLA 필라멘트의 물성을 제공하고 있지는 않는 것 같음.
일단은 3D-printing된 PLA로 tensile test 진행한 [논문](https://www.mdpi.com/2073-4360/12/1/250)의 stress-strain 그래프에서 점 찍어서 좌표 추출해서 ``PLA stress-strain.xlsx``에 저장했고, 이를 바탕으로 시뮬레이션을 진행하고 있음.

![Stress-Strain of PLA](https://www.mdpi.com/polymers/polymers-12-00250/article_deploy/html/images/polymers-12-00250-g005-550.jpg)
