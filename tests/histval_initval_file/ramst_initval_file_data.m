x = vertcat([ 1; 1.2 ], repmat(1, 200, 1));
k = repmat(13, 202, 1);
c = repmat(1.5, 202, 1);
save('ramst_initval_file_data_col_vec_mat.mat','c','k','x');

if ispc()
    xlswrite('ramst_initval_file_excel',[x k c],1,'A2');
    xlswrite('ramst_initval_file_excel',{'x' 'k' 'c'},1,'A1');
end

c=c';
k=k';
x=x';
save('ramst_initval_file_data_row_vec_mat.mat','c','k','x');
