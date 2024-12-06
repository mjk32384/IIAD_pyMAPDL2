## 파일 설명
``functions.py``: tensile test 진행 코드, CAD 생성 코드 있는 파일.  
``test_*.py``: ``functions.py``에 적절한 파라미터를 넘겨주어서, 원하는 시뮬레이션을 하는 파일.  
``result_visulization.ipynb``: 시뮬레이션 결과를 시각화한 파일. 그래프 등등 다 여기 있음. 근데 mapdl 인스턴스를 이용해서 시각화한거는 다른 컴퓨터에서는 안 보일 거임.  
``tensile_*.ipynb``: 이것저것 시도해본 파일. 무시해도 됨.

``result/``: 시뮬레이션 결과가 저장되는 폴더. 용량 문제로 git에 업로드 안 하고 있음.  
``Examples/``: pyMAPDL 예제 있는 폴더. 무시해도 됨.

## 시뮬레이션 방법
Damage/Fracture Model 혹은 failure criterion에 대해서 조사하고, 무엇이 3D 프린팅된 PLA에 적합한지 결정해야 할 듯.  
현재는 Multilinear Isotropic hardening model (MISO Model)을 통해 non-linear 부분을 시뮬레이션 하고, Maximum Normal Stress Theory를 이용해서 fracture 여부를 판단하고 있음.

[신도리코](https://sindoh4u.com/goods/view?no=46&srsltid=AfmBOopVEhJF0xLJstQIMS8SeHpXZjQAaL1keqCtixq9tZ3yZEf-YIp9)에서는 PLA 필라멘트의 물성을 제공하고 있지는 않는 것 같음.
일단은 3D-printing된 PLA로 tensile test 진행한 [논문](https://www.mdpi.com/2073-4360/12/1/250)의 stress-strain 그래프에서 점 찍어서 좌표 추출해서 ``PLA stress-strain.xlsx``에 저장했음.

![Stress-Strain of PLA](https://www.mdpi.com/polymers/polymers-12-00250/article_deploy/html/images/polymers-12-00250-g005-550.jpg)
