% Liet ke tat ca cac so nguyen to tu 0 - 1000.
% Y tuong: Viet 1 ham kiem tra 1 so co phai nguyen to hay khong.
% Sau do chay vong lap, neu dung la so nguyen to thi cho vao 1 array
% Cuoi cung in array do ra.

function [value] = kiemtra_ngto(n)

k = ceil(sqrt(n));

value = 1; % Ban dau coi nhu n la so nguyen to

if n >= 3
    for i = 2:k
        if mod(n,i) == 0
            value = 0;
            break
        end
    end
end
