load exampleData
k = data(1:2:end) + 1i*(data(2:2:end));
k = reshape(k, [128 15 128]);
k = permute(k, [1 3 2]);
imagesc(abs(k(:,:,1)))
imagesc(log(abs(k(:,:,1))))
IM = abs(  fftshift(  ifft2(k)));
imagesc(IM(:,:,1));
colormap gray
axis image
imagesc(IM(:,:,10));
colormap gray
axis image