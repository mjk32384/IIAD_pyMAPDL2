``tensile_test.ipynb``: tensile test를 비롯한 이것저것 테스트해보고 있는 파일  
``result``: 시뮬레이션 결과가 저장되는 폴더. 나중에는 용량 문제로 git에 저장 안 할 수도 있을 듯.

Damage/Fracture Model 혹은 failure criterion에 대해서 조사하고, 무엇이 3D 프린팅된 PLA에 적합한지 결정해야 할 듯.  
현재는 Multilinear Isotropic hardening model (MISO Model)을 통해 non-linear 부분을 시뮬레이션 하고, Maximum Normal Stress Theory를 이용해서 fracture 여부를 판단하고 있음.

[신도리코](https://sindoh4u.com/goods/view?no=46&srsltid=AfmBOopVEhJF0xLJstQIMS8SeHpXZjQAaL1keqCtixq9tZ3yZEf-YIp9)에서는 PLA 필라멘트의 물성을 제공하고 있지는 않는 것 같음.
일단은 3D-printing된 PLA로 tensile test 진행한 [논문](https://www.mdpi.com/2073-4360/12/1/250)의 stress-strain 그래프에서 점 찍어서 좌표 추출해서 물성값 정했음.  
![Stress-Strain of PLA](https://www.mdpi.com/polymers/polymers-12-00250/article_deploy/html/images/polymers-12-00250-g005-550.jpg)
