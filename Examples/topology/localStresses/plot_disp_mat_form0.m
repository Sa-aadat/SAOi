ncol=15*2;
nrow= 5*2;

load plotent.out

x_mat=zeros(nrow,ncol);

for icol=1:ncol
  pent_rows=(1:nrow) + (icol-1)*nrow;
  x_mat(:,icol)=1-plotent(pent_rows,1);
end

%x_mat=[fliplr(x_mat) x_mat];

colormap('gray')
imagesc( x_mat )
axis equal
axis off
