VBControl - это бибилотека на C++, реализующая основные операции для проектирования систем управления. Также имеет python-обертку. <br>

## Возможности:
 - преобразование непрерывных линейных систем в дискретные
 - ПИД-регулятор (c функцией устарения ошибки)
 - П-регулятор полного состояния системы
 - Модель балансирующего робота (линейная+нелинейная+якобиан)
 - Модель DC-мотора (линейная + нелинейная)
 - Модель PMSM-мотора с синусоидальным противо-ЭДС и учтенным зубцовым моментом (нелинейная)
 - FOC для PMSM
 - Решение дискретной линейно-квадратичной задачи оптимального управления
 - Low Pass Filter
 - Linear Kalman Filter
 - Extended Kalman Filter
 - Unscented Kalman Filter
 - Разные математические действия (clip, sign, friction compensation)

## Зависимости
 - pybind11
 - matploplusplus
 - Eigen3

## Установка
```
cd ~
git clone https://github.com/voltbro/vbcontrol.git
cd vbcontrol
mkdir build
cd build
cmake ..
make
```