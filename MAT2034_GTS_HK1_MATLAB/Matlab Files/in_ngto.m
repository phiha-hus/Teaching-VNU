% In cac so nguyen to tu 1 den 1000000 va dem thoi gian
tic
S = [];
for n = 2:1e+6
    if kiemtra_ngto(n) == 1
        S = [S n];
    end
end
toc

disp(strcat('Thoi gian de tim tat ca cac so ngto tu 1 den 10000 la: ',num2str(toc)))

