function y         	= dpm_interpn(varargin)
switch((nargin-1)/2)
    case 1
        xx1 = varargin{1};
        YY  = varargin{2};
        A1  = varargin{3};

        Ars = [reshape(A1, [numel(A1) ,1])];
        xx{1}  = reshape(xx1,[numel(xx1),1]);

        h = max([max(diff(xx{1}))],eps);

        Ars(Ars(:,1)<min(xx{1}),1) = min(xx{1});

        Ars(Ars(:,1)>max(xx{1}),1) = max(xx{1});

        ind(:,1,1)  = 1 + floor(round((Ars(:,1)-xx{1}(1))/h(1)*1e8)*1e-8);

        ind(:,1,2) = 1 +  ceil(round((Ars(:,1)-xx{1}(1))/h(1)*1e8)*1e-8);

        yy(:,1,1)    = YY(dpm_sub2ind(size(YY),ind(:,1,1)));
        yy(:,2,1)    = YY(dpm_sub2ind(size(YY),ind(:,1,2)));

        da(:,1) = (Ars(:,1)-xx{1}(ind(:,1,1)))/h(1);

        yi = yy(:,1,1) + da(:,1).*(yy(:,2,1) - yy(:,1,1));

        y = reshape(yi,size(A1));

    case 2
        xx2 = varargin{1};
        xx1 = varargin{2};
        YY  = varargin{3};
        A2  = varargin{4};
        A1  = varargin{5};

        Ars = [reshape(A1, [numel(A1) ,1]) reshape(A2, [numel(A2) ,1])];
        xx{1}  = reshape(xx1,[numel(xx1),1]);
        xx{2}  = reshape(xx2,[numel(xx2),1]);

        h = max([max(diff(xx{1})) max(diff(xx{2}))],eps);

        Ars(Ars(:,1)<min(xx{1}),1) = min(xx{1});
        Ars(Ars(:,2)<min(xx{2}),2) = min(xx{2});

        Ars(Ars(:,1)>max(xx{1}),1) = max(xx{1});
        Ars(Ars(:,2)>max(xx{2}),2) = max(xx{2});

        ind(:,1,1)  = 1 + floor(round((Ars(:,1)-xx{1}(1))/h(1)*1e8)*1e-8);
        ind(:,2,1)  = 1 + floor(round((Ars(:,2)-xx{2}(1))/h(2)*1e8)*1e-8);
        ind(:,1,2) = 1 +  ceil(round((Ars(:,1)-xx{1}(1))/h(1)*1e8)*1e-8);
        ind(:,2,2) = 1 +  ceil(round((Ars(:,2)-xx{2}(1))/h(2)*1e8)*1e-8);

        yy(:,1,1)    = YY(dpm_sub2ind(size(YY),ind(:,1,1),ind(:,2,1)));
        yy(:,2,1)    = YY(dpm_sub2ind(size(YY),ind(:,1,2),ind(:,2,1)));
        yy(:,1,2)    = YY(dpm_sub2ind(size(YY),ind(:,1,1),ind(:,2,2)));
        yy(:,2,2)    = YY(dpm_sub2ind(size(YY),ind(:,1,2),ind(:,2,2)));

        da(:,1) = (Ars(:,1)-xx{1}(ind(:,1,1)))/h(1);
        da(:,2) = (Ars(:,2)-xx{2}(ind(:,2,1)))/h(2);

        v1(:,1) = yy(:,1,1) + da(:,1).*(yy(:,2,1) - yy(:,1,1));
        v1(:,2) = yy(:,1,2) + da(:,1).*(yy(:,2,2) - yy(:,1,2));

        yi = da(:,2).*(v1(:,2) - v1(:,1)) + v1(:,1);

        y = reshape(yi,size(A1));
    case 3
        xx3 = varargin{1};
        xx2 = varargin{2};
        xx1 = varargin{3};
        YY  = varargin{4};
        A3  = varargin{5};
        A2  = varargin{6};
        A1  = varargin{7};

        Ars = [reshape(A1, [numel(A1) ,1]) reshape(A2, [numel(A2) ,1]) reshape(A3, [numel(A3) ,1])];
        xx{1}  = reshape(xx1,[numel(xx1),1]);
        xx{2}  = reshape(xx2,[numel(xx2),1]);
        xx{3}  = reshape(xx3,[numel(xx3),1]);

        h = max([max(diff(xx{1})) max(diff(xx{2})) max(diff(xx{3}))],eps);

        Ars(Ars(:,1)<min(xx{1}),1) = min(xx{1});
        Ars(Ars(:,2)<min(xx{2}),2) = min(xx{2});
        Ars(Ars(:,3)<min(xx{3}),3) = min(xx{3});

        Ars(Ars(:,1)>max(xx{1}),1) = max(xx{1});
        Ars(Ars(:,2)>max(xx{2}),2) = max(xx{2});
        Ars(Ars(:,3)>max(xx{3}),3) = max(xx{3});

        ind(:,1,1)  = 1 + floor(round((Ars(:,1)-xx{1}(1))/h(1)*1e8)*1e-8);
        ind(:,2,1)  = 1 + floor(round((Ars(:,2)-xx{2}(1))/h(2)*1e8)*1e-8);
        ind(:,3,1)  = 1 + floor(round((Ars(:,3)-xx{3}(1))/h(3)*1e8)*1e-8);

        ind(:,1,2) = 1 +  ceil(round((Ars(:,1)-xx{1}(1))/h(1)*1e8)*1e-8);
        ind(:,2,2) = 1 +  ceil(round((Ars(:,2)-xx{2}(1))/h(2)*1e8)*1e-8);
        ind(:,3,2) = 1 +  ceil(round((Ars(:,3)-xx{3}(1))/h(3)*1e8)*1e-8);

        yy(:,1,1,1)    = YY(dpm_sub2ind(size(YY),ind(:,1,1),ind(:,2,1),ind(:,3,1)));
        yy(:,2,1,1)    = YY(dpm_sub2ind(size(YY),ind(:,1,2),ind(:,2,1),ind(:,3,1)));
        yy(:,1,2,1)    = YY(dpm_sub2ind(size(YY),ind(:,1,1),ind(:,2,2),ind(:,3,1)));
        yy(:,2,2,1)    = YY(dpm_sub2ind(size(YY),ind(:,1,2),ind(:,2,2),ind(:,3,1)));
        yy(:,1,1,2)    = YY(dpm_sub2ind(size(YY),ind(:,1,1),ind(:,2,1),ind(:,3,2)));
        yy(:,2,1,2)    = YY(dpm_sub2ind(size(YY),ind(:,1,2),ind(:,2,1),ind(:,3,2)));
        yy(:,1,2,2)    = YY(dpm_sub2ind(size(YY),ind(:,1,1),ind(:,2,2),ind(:,3,2)));
        yy(:,2,2,2)    = YY(dpm_sub2ind(size(YY),ind(:,1,2),ind(:,2,2),ind(:,3,2)));

        da(:,1) = (Ars(:,1)-xx{1}(ind(:,1,1)))/h(1);
        da(:,2) = (Ars(:,2)-xx{2}(ind(:,2,1)))/h(2);
        da(:,3) = (Ars(:,3)-xx{3}(ind(:,3,1)))/h(3);

        
        v2(:,1,1) = yy(:,1,1,1) + da(:,3).*(yy(:,1,1,2) - yy(:,1,1,1));
        v2(:,2,1) = yy(:,2,1,1) + da(:,3).*(yy(:,2,1,2) - yy(:,2,1,1));
        v2(:,1,2) = yy(:,1,2,1) + da(:,3).*(yy(:,1,2,2) - yy(:,1,2,1));
        v2(:,2,2) = yy(:,2,2,1) + da(:,3).*(yy(:,2,2,2) - yy(:,2,2,1)); 
        
        v1(:,1) = v2(:,1,1) + da(:,2).*(v2(:,1,2) - v2(:,1,1));
        v1(:,2) = v2(:,2,1) + da(:,2).*(v2(:,2,2) - v2(:,2,1));

        yi = da(:,1).*(v1(:,2) - v1(:,1)) + v1(:,1);        
        
        
        
        
        
        y = reshape(yi,size(A1));
    case 4
        xx4 = varargin{1};
        xx3 = varargin{2};
        xx2 = varargin{3};
        xx1 = varargin{4};
        YY  = varargin{5};
        A4  = varargin{6};
        A3  = varargin{7};
        A2  = varargin{8};
        A1  = varargin{9};

        Ars = [reshape(A1, [numel(A1) ,1]) reshape(A2, [numel(A2) ,1]) reshape(A3, [numel(A3) ,1]) reshape(A4, [numel(A4) ,1])];
        xx{1}  = reshape(xx1,[numel(xx1),1]);
        xx{2}  = reshape(xx2,[numel(xx2),1]);
        xx{3}  = reshape(xx3,[numel(xx3),1]);
        xx{4}  = reshape(xx4,[numel(xx4),1]);

        h = max([max(diff(xx{1})) max(diff(xx{2})) max(diff(xx{3})) max(diff(xx{4}))],eps);

        Ars(Ars(:,1)<min(xx{1}),1) = min(xx{1});
        Ars(Ars(:,2)<min(xx{2}),2) = min(xx{2});
        Ars(Ars(:,3)<min(xx{3}),3) = min(xx{3});
        Ars(Ars(:,4)<min(xx{4}),4) = min(xx{4});

        Ars(Ars(:,1)>max(xx{1}),1) = max(xx{1});
        Ars(Ars(:,2)>max(xx{2}),2) = max(xx{2});
        Ars(Ars(:,3)>max(xx{3}),3) = max(xx{3});
        Ars(Ars(:,4)>max(xx{4}),4) = max(xx{4});

        ind(:,1,1)  = 1 + floor(round((Ars(:,1)-xx{1}(1))/h(1)*1e8)*1e-8);
        ind(:,2,1)  = 1 + floor(round((Ars(:,2)-xx{2}(1))/h(2)*1e8)*1e-8);
        ind(:,3,1)  = 1 + floor(round((Ars(:,3)-xx{3}(1))/h(3)*1e8)*1e-8);
        ind(:,4,1)  = 1 + floor(round((Ars(:,4)-xx{4}(1))/h(4)*1e8)*1e-8);

        ind(:,1,2) = 1 +  ceil(round((Ars(:,1)-xx{1}(1))/h(1)*1e8)*1e-8);
        ind(:,2,2) = 1 +  ceil(round((Ars(:,2)-xx{2}(1))/h(2)*1e8)*1e-8);
        ind(:,3,2) = 1 +  ceil(round((Ars(:,3)-xx{3}(1))/h(3)*1e8)*1e-8);
        ind(:,4,2) = 1 +  ceil(round((Ars(:,4)-xx{4}(1))/h(4)*1e8)*1e-8);

        yy(:,1,1,1,1)    = YY(dpm_sub2ind(size(YY),ind(:,1,1),ind(:,2,1),ind(:,3,1),ind(:,4,1)));
        yy(:,2,1,1,1)    = YY(dpm_sub2ind(size(YY),ind(:,1,2),ind(:,2,1),ind(:,3,1),ind(:,4,1)));
        yy(:,1,2,1,1)    = YY(dpm_sub2ind(size(YY),ind(:,1,1),ind(:,2,2),ind(:,3,1),ind(:,4,1)));
        yy(:,2,2,1,1)    = YY(dpm_sub2ind(size(YY),ind(:,1,2),ind(:,2,2),ind(:,3,1),ind(:,4,1)));
        yy(:,1,1,2,1)    = YY(dpm_sub2ind(size(YY),ind(:,1,1),ind(:,2,1),ind(:,3,2),ind(:,4,1)));
        yy(:,2,1,2,1)    = YY(dpm_sub2ind(size(YY),ind(:,1,2),ind(:,2,1),ind(:,3,2),ind(:,4,1)));
        yy(:,1,2,2,1)    = YY(dpm_sub2ind(size(YY),ind(:,1,1),ind(:,2,2),ind(:,3,2),ind(:,4,1)));
        yy(:,2,2,2,1)    = YY(dpm_sub2ind(size(YY),ind(:,1,2),ind(:,2,2),ind(:,3,2),ind(:,4,1)));
        yy(:,1,1,1,2)    = YY(dpm_sub2ind(size(YY),ind(:,1,1),ind(:,2,1),ind(:,3,1),ind(:,4,2)));
        yy(:,2,1,1,2)    = YY(dpm_sub2ind(size(YY),ind(:,1,2),ind(:,2,1),ind(:,3,1),ind(:,4,2)));
        yy(:,1,2,1,2)    = YY(dpm_sub2ind(size(YY),ind(:,1,1),ind(:,2,2),ind(:,3,1),ind(:,4,2)));
        yy(:,2,2,1,2)    = YY(dpm_sub2ind(size(YY),ind(:,1,2),ind(:,2,2),ind(:,3,1),ind(:,4,2)));
        yy(:,1,1,2,2)    = YY(dpm_sub2ind(size(YY),ind(:,1,1),ind(:,2,1),ind(:,3,2),ind(:,4,2)));
        yy(:,2,1,2,2)    = YY(dpm_sub2ind(size(YY),ind(:,1,2),ind(:,2,1),ind(:,3,2),ind(:,4,2)));
        yy(:,1,2,2,2)    = YY(dpm_sub2ind(size(YY),ind(:,1,1),ind(:,2,2),ind(:,3,2),ind(:,4,2)));
        yy(:,2,2,2,2)    = YY(dpm_sub2ind(size(YY),ind(:,1,2),ind(:,2,2),ind(:,3,2),ind(:,4,2)));

        da(:,1) = (Ars(:,1)-xx{1}(ind(:,1,1)))/h(1);
        da(:,2) = (Ars(:,2)-xx{2}(ind(:,2,1)))/h(2);
        da(:,3) = (Ars(:,3)-xx{3}(ind(:,3,1)))/h(3);
        da(:,4) = (Ars(:,4)-xx{4}(ind(:,4,1)))/h(4);

        
        v3(:,1,1,1) = yy(:,1,1,1,1) + da(:,4).*(yy(:,1,1,1,2) - yy(:,1,1,1,1));
        v3(:,2,1,1) = yy(:,2,1,1,1) + da(:,4).*(yy(:,2,1,1,2) - yy(:,2,1,1,1));
        v3(:,1,2,1) = yy(:,1,2,1,1) + da(:,4).*(yy(:,1,2,1,2) - yy(:,1,2,1,1));
        v3(:,2,2,1) = yy(:,2,2,1,1) + da(:,4).*(yy(:,2,2,1,2) - yy(:,2,2,1,1));
        v3(:,1,1,2) = yy(:,1,1,2,1) + da(:,4).*(yy(:,1,1,2,2) - yy(:,1,1,2,1));
        v3(:,2,1,2) = yy(:,2,1,2,1) + da(:,4).*(yy(:,2,1,2,2) - yy(:,2,1,2,1));
        v3(:,1,2,2) = yy(:,1,2,2,1) + da(:,4).*(yy(:,1,2,2,2) - yy(:,1,2,2,1));
        v3(:,2,2,2) = yy(:,2,2,2,1) + da(:,4).*(yy(:,2,2,2,2) - yy(:,2,2,2,1));
        
        v2(:,1,1) = v3(:,1,1,1) + da(:,3).*(v3(:,1,1,2) - v3(:,1,1,1));
        v2(:,2,1) = v3(:,2,1,1) + da(:,3).*(v3(:,2,1,2) - v3(:,2,1,1));
        v2(:,1,2) = v3(:,1,2,1) + da(:,3).*(v3(:,1,2,2) - v3(:,1,2,1));
        v2(:,2,2) = v3(:,2,2,1) + da(:,3).*(v3(:,2,2,2) - v3(:,2,2,1)); 
        
        v1(:,1) = v2(:,1,1) + da(:,2).*(v2(:,1,2) - v2(:,1,1));
        v1(:,2) = v2(:,2,1) + da(:,2).*(v2(:,2,2) - v2(:,2,1));

        yi = da(:,1).*(v1(:,2) - v1(:,1)) + v1(:,1);          
        
        
        
        y = reshape(yi,size(A1));
    case 5
        xx5 = varargin{1};
        xx4 = varargin{2};
        xx3 = varargin{3};
        xx2 = varargin{4};
        xx1 = varargin{5};
        YY  = varargin{6};
        A5  = varargin{7};
        A4  = varargin{8};
        A3  = varargin{9};
        A2  = varargin{10};
        A1  = varargin{11};

        Ars = [reshape(A1, [numel(A1) ,1]) reshape(A2, [numel(A2) ,1]) reshape(A3, [numel(A3) ,1]) reshape(A4, [numel(A4) ,1]) reshape(A5, [numel(A5) ,1])];
        xx{1}  = reshape(xx1,[numel(xx1),1]);
        xx{2}  = reshape(xx2,[numel(xx2),1]);
        xx{3}  = reshape(xx3,[numel(xx3),1]);
        xx{4}  = reshape(xx4,[numel(xx4),1]);
        xx{5}  = reshape(xx5,[numel(xx5),1]);

        h = max([max(diff(xx{1})) max(diff(xx{2})) max(diff(xx{3})) max(diff(xx{4})) max(diff(xx{5}))],eps);

        Ars(Ars(:,1)<min(xx{1}),1) = min(xx{1});
        Ars(Ars(:,2)<min(xx{2}),2) = min(xx{2});
        Ars(Ars(:,3)<min(xx{3}),3) = min(xx{3});
        Ars(Ars(:,4)<min(xx{4}),4) = min(xx{4});
        Ars(Ars(:,5)<min(xx{5}),5) = min(xx{5});

        Ars(Ars(:,1)>max(xx{1}),1) = max(xx{1});
        Ars(Ars(:,2)>max(xx{2}),2) = max(xx{2});
        Ars(Ars(:,3)>max(xx{3}),3) = max(xx{3});
        Ars(Ars(:,4)>max(xx{4}),4) = max(xx{4});
        Ars(Ars(:,5)>max(xx{5}),5) = max(xx{5});

        ind(:,1,1)  = 1 + floor(round((Ars(:,1)-xx{1}(1))/h(1)*1e8)*1e-8);
        ind(:,2,1)  = 1 + floor(round((Ars(:,2)-xx{2}(1))/h(2)*1e8)*1e-8);
        ind(:,3,1)  = 1 + floor(round((Ars(:,3)-xx{3}(1))/h(3)*1e8)*1e-8);
        ind(:,4,1)  = 1 + floor(round((Ars(:,4)-xx{4}(1))/h(4)*1e8)*1e-8);
        ind(:,5,1)  = 1 + floor(round((Ars(:,5)-xx{5}(1))/h(5)*1e8)*1e-8);

        ind(:,1,2) = 1 +  ceil(round((Ars(:,1)-xx{1}(1))/h(1)*1e8)*1e-8);
        ind(:,2,2) = 1 +  ceil(round((Ars(:,2)-xx{2}(1))/h(2)*1e8)*1e-8);
        ind(:,3,2) = 1 +  ceil(round((Ars(:,3)-xx{3}(1))/h(3)*1e8)*1e-8);
        ind(:,4,2) = 1 +  ceil(round((Ars(:,4)-xx{4}(1))/h(4)*1e8)*1e-8);
        ind(:,5,2) = 1 +  ceil(round((Ars(:,5)-xx{5}(1))/h(5)*1e8)*1e-8);

        yy(:,1,1,1,1,1)    = YY(dpm_sub2ind(size(YY),ind(:,1,1),ind(:,2,1),ind(:,3,1),ind(:,4,1),ind(:,5,1)));
        yy(:,2,1,1,1,1)    = YY(dpm_sub2ind(size(YY),ind(:,1,2),ind(:,2,1),ind(:,3,1),ind(:,4,1),ind(:,5,1)));
        yy(:,1,2,1,1,1)    = YY(dpm_sub2ind(size(YY),ind(:,1,1),ind(:,2,2),ind(:,3,1),ind(:,4,1),ind(:,5,1)));
        yy(:,2,2,1,1,1)    = YY(dpm_sub2ind(size(YY),ind(:,1,2),ind(:,2,2),ind(:,3,1),ind(:,4,1),ind(:,5,1)));
        yy(:,1,1,2,1,1)    = YY(dpm_sub2ind(size(YY),ind(:,1,1),ind(:,2,1),ind(:,3,2),ind(:,4,1),ind(:,5,1)));
        yy(:,2,1,2,1,1)    = YY(dpm_sub2ind(size(YY),ind(:,1,2),ind(:,2,1),ind(:,3,2),ind(:,4,1),ind(:,5,1)));
        yy(:,1,2,2,1,1)    = YY(dpm_sub2ind(size(YY),ind(:,1,1),ind(:,2,2),ind(:,3,2),ind(:,4,1),ind(:,5,1)));
        yy(:,2,2,2,1,1)    = YY(dpm_sub2ind(size(YY),ind(:,1,2),ind(:,2,2),ind(:,3,2),ind(:,4,1),ind(:,5,1)));
        yy(:,1,1,1,2,1)    = YY(dpm_sub2ind(size(YY),ind(:,1,1),ind(:,2,1),ind(:,3,1),ind(:,4,2),ind(:,5,1)));
        yy(:,2,1,1,2,1)    = YY(dpm_sub2ind(size(YY),ind(:,1,2),ind(:,2,1),ind(:,3,1),ind(:,4,2),ind(:,5,1)));
        yy(:,1,2,1,2,1)    = YY(dpm_sub2ind(size(YY),ind(:,1,1),ind(:,2,2),ind(:,3,1),ind(:,4,2),ind(:,5,1)));
        yy(:,2,2,1,2,1)    = YY(dpm_sub2ind(size(YY),ind(:,1,2),ind(:,2,2),ind(:,3,1),ind(:,4,2),ind(:,5,1)));
        yy(:,1,1,2,2,1)    = YY(dpm_sub2ind(size(YY),ind(:,1,1),ind(:,2,1),ind(:,3,2),ind(:,4,2),ind(:,5,1)));
        yy(:,2,1,2,2,1)    = YY(dpm_sub2ind(size(YY),ind(:,1,2),ind(:,2,1),ind(:,3,2),ind(:,4,2),ind(:,5,1)));
        yy(:,1,2,2,2,1)    = YY(dpm_sub2ind(size(YY),ind(:,1,1),ind(:,2,2),ind(:,3,2),ind(:,4,2),ind(:,5,1)));
        yy(:,2,2,2,2,1)    = YY(dpm_sub2ind(size(YY),ind(:,1,2),ind(:,2,2),ind(:,3,2),ind(:,4,2),ind(:,5,1)));
        yy(:,1,1,1,1,2)    = YY(dpm_sub2ind(size(YY),ind(:,1,1),ind(:,2,1),ind(:,3,1),ind(:,4,1),ind(:,5,2)));
        yy(:,2,1,1,1,2)    = YY(dpm_sub2ind(size(YY),ind(:,1,2),ind(:,2,1),ind(:,3,1),ind(:,4,1),ind(:,5,2)));
        yy(:,1,2,1,1,2)    = YY(dpm_sub2ind(size(YY),ind(:,1,1),ind(:,2,2),ind(:,3,1),ind(:,4,1),ind(:,5,2)));
        yy(:,2,2,1,1,2)    = YY(dpm_sub2ind(size(YY),ind(:,1,2),ind(:,2,2),ind(:,3,1),ind(:,4,1),ind(:,5,2)));
        yy(:,1,1,2,1,2)    = YY(dpm_sub2ind(size(YY),ind(:,1,1),ind(:,2,1),ind(:,3,2),ind(:,4,1),ind(:,5,2)));
        yy(:,2,1,2,1,2)    = YY(dpm_sub2ind(size(YY),ind(:,1,2),ind(:,2,1),ind(:,3,2),ind(:,4,1),ind(:,5,2)));
        yy(:,1,2,2,1,2)    = YY(dpm_sub2ind(size(YY),ind(:,1,1),ind(:,2,2),ind(:,3,2),ind(:,4,1),ind(:,5,2)));
        yy(:,2,2,2,1,2)    = YY(dpm_sub2ind(size(YY),ind(:,1,2),ind(:,2,2),ind(:,3,2),ind(:,4,1),ind(:,5,2)));
        yy(:,1,1,1,2,2)    = YY(dpm_sub2ind(size(YY),ind(:,1,1),ind(:,2,1),ind(:,3,1),ind(:,4,2),ind(:,5,2)));
        yy(:,2,1,1,2,2)    = YY(dpm_sub2ind(size(YY),ind(:,1,2),ind(:,2,1),ind(:,3,1),ind(:,4,2),ind(:,5,2)));
        yy(:,1,2,1,2,2)    = YY(dpm_sub2ind(size(YY),ind(:,1,1),ind(:,2,2),ind(:,3,1),ind(:,4,2),ind(:,5,2)));
        yy(:,2,2,1,2,2)    = YY(dpm_sub2ind(size(YY),ind(:,1,2),ind(:,2,2),ind(:,3,1),ind(:,4,2),ind(:,5,2)));
        yy(:,1,1,2,2,2)    = YY(dpm_sub2ind(size(YY),ind(:,1,1),ind(:,2,1),ind(:,3,2),ind(:,4,2),ind(:,5,2)));
        yy(:,2,1,2,2,2)    = YY(dpm_sub2ind(size(YY),ind(:,1,2),ind(:,2,1),ind(:,3,2),ind(:,4,2),ind(:,5,2)));
        yy(:,1,2,2,2,2)    = YY(dpm_sub2ind(size(YY),ind(:,1,1),ind(:,2,2),ind(:,3,2),ind(:,4,2),ind(:,5,2)));
        yy(:,2,2,2,2,2)    = YY(dpm_sub2ind(size(YY),ind(:,1,2),ind(:,2,2),ind(:,3,2),ind(:,4,2),ind(:,5,2)));

        da(:,1) = (Ars(:,1)-xx{1}(ind(:,1,1)))/h(1);
        da(:,2) = (Ars(:,2)-xx{2}(ind(:,2,1)))/h(2);
        da(:,3) = (Ars(:,3)-xx{3}(ind(:,3,1)))/h(3);
        da(:,4) = (Ars(:,4)-xx{4}(ind(:,4,1)))/h(4);
        da(:,5) = (Ars(:,5)-xx{5}(ind(:,5,1)))/h(5);

        v4(:,1,1,1,1) = yy(:,1,1,1,1,1) + da(:,5).*(yy(:,1,1,1,1,2) - yy(:,1,1,1,1,1));
        v4(:,2,1,1,1) = yy(:,2,1,1,1,1) + da(:,5).*(yy(:,2,1,1,1,2) - yy(:,2,1,1,1,1));
        v4(:,1,2,1,1) = yy(:,1,2,1,1,1) + da(:,5).*(yy(:,1,2,1,1,2) - yy(:,1,2,1,1,1));
        v4(:,2,2,1,1) = yy(:,2,2,1,1,1) + da(:,5).*(yy(:,2,2,1,1,2) - yy(:,2,2,1,1,1));
        v4(:,1,1,2,1) = yy(:,1,1,2,1,1) + da(:,5).*(yy(:,1,1,2,1,2) - yy(:,1,1,2,1,1));
        v4(:,2,1,2,1) = yy(:,2,1,2,1,1) + da(:,5).*(yy(:,2,1,2,1,2) - yy(:,2,1,2,1,1));
        v4(:,1,2,2,1) = yy(:,1,2,2,1,1) + da(:,5).*(yy(:,1,2,2,1,2) - yy(:,1,2,2,1,1));
        v4(:,2,2,2,1) = yy(:,2,2,2,1,1) + da(:,5).*(yy(:,2,2,2,1,2) - yy(:,2,2,2,1,1));
        v4(:,1,1,1,2) = yy(:,1,1,1,2,1) + da(:,5).*(yy(:,1,1,1,2,2) - yy(:,1,1,1,2,1));
        v4(:,2,1,1,2) = yy(:,2,1,1,2,1) + da(:,5).*(yy(:,2,1,1,2,2) - yy(:,2,1,1,2,1));
        v4(:,1,2,1,2) = yy(:,1,2,1,2,1) + da(:,5).*(yy(:,1,2,1,2,2) - yy(:,1,2,1,2,1));
        v4(:,2,2,1,2) = yy(:,2,2,1,2,1) + da(:,5).*(yy(:,2,2,1,2,2) - yy(:,2,2,1,2,1));
        v4(:,1,1,2,2) = yy(:,1,1,2,2,1) + da(:,5).*(yy(:,1,1,2,2,2) - yy(:,1,1,2,2,1));
        v4(:,2,1,2,2) = yy(:,2,1,2,2,1) + da(:,5).*(yy(:,2,1,2,2,2) - yy(:,2,1,2,2,1));
        v4(:,1,2,2,2) = yy(:,1,2,2,2,1) + da(:,5).*(yy(:,1,2,2,2,2) - yy(:,1,2,2,2,1));
        v4(:,2,2,2,2) = yy(:,2,2,2,2,1) + da(:,5).*(yy(:,2,2,2,2,2) - yy(:,2,2,2,2,1));


        v3(:,1,1,1) = v4(:,1,1,1,1) + da(:,4).*(v4(:,1,1,1,2) - v4(:,1,1,1,1));
        v3(:,2,1,1) = v4(:,2,1,1,1) + da(:,4).*(v4(:,2,1,1,2) - v4(:,2,1,1,1));
        v3(:,1,2,1) = v4(:,1,2,1,1) + da(:,4).*(v4(:,1,2,1,2) - v4(:,1,2,1,1));
        v3(:,2,2,1) = v4(:,2,2,1,1) + da(:,4).*(v4(:,2,2,1,2) - v4(:,2,2,1,1));
        v3(:,1,1,2) = v4(:,1,1,2,1) + da(:,4).*(v4(:,1,1,2,2) - v4(:,1,1,2,1));
        v3(:,2,1,2) = v4(:,2,1,2,1) + da(:,4).*(v4(:,2,1,2,2) - v4(:,2,1,2,1));
        v3(:,1,2,2) = v4(:,1,2,2,1) + da(:,4).*(v4(:,1,2,2,2) - v4(:,1,2,2,1));
        v3(:,2,2,2) = v4(:,2,2,2,1) + da(:,4).*(v4(:,2,2,2,2) - v4(:,2,2,2,1));
        
        v2(:,1,1) = v3(:,1,1,1) + da(:,3).*(v3(:,1,1,2) - v3(:,1,1,1));
        v2(:,2,1) = v3(:,2,1,1) + da(:,3).*(v3(:,2,1,2) - v3(:,2,1,1));
        v2(:,1,2) = v3(:,1,2,1) + da(:,3).*(v3(:,1,2,2) - v3(:,1,2,1));
        v2(:,2,2) = v3(:,2,2,1) + da(:,3).*(v3(:,2,2,2) - v3(:,2,2,1)); 
        
        v1(:,1) = v2(:,1,1) + da(:,2).*(v2(:,1,2) - v2(:,1,1));
        v1(:,2) = v2(:,2,1) + da(:,2).*(v2(:,2,2) - v2(:,2,1));

        yi = da(:,1).*(v1(:,2) - v1(:,1)) + v1(:,1);            
        
        y = reshape(yi,size(A1));
    case 6
        xx6 = varargin{1};
        xx5 = varargin{2};
        xx4 = varargin{3};
        xx3 = varargin{4};
        xx2 = varargin{5};
        xx1 = varargin{6};
        YY  = varargin{7};
        A6  = varargin{8};
        A5  = varargin{9};
        A4  = varargin{10};
        A3  = varargin{11};
        A2  = varargin{12};
        A1  = varargin{13};

        Ars = [reshape(A1, [numel(A1) ,1]) reshape(A2, [numel(A2) ,1]) reshape(A3, [numel(A3) ,1]) reshape(A4, [numel(A4) ,1]) reshape(A5, [numel(A5) ,1]) reshape(A6, [numel(A6) ,1])];
        xx{1}  = reshape(xx1,[numel(xx1),1]);
        xx{2}  = reshape(xx2,[numel(xx2),1]);
        xx{3}  = reshape(xx3,[numel(xx3),1]);
        xx{4}  = reshape(xx4,[numel(xx4),1]);
        xx{5}  = reshape(xx5,[numel(xx5),1]);
        xx{6}  = reshape(xx6,[numel(xx6),1]);

        h = max([max(diff(xx{1})) max(diff(xx{2})) max(diff(xx{3})) max(diff(xx{4})) max(diff(xx{5})) max(diff(xx{6}))],eps);

        Ars(Ars(:,1)<min(xx{1}),1) = min(xx{1});
        Ars(Ars(:,2)<min(xx{2}),2) = min(xx{2});
        Ars(Ars(:,3)<min(xx{3}),3) = min(xx{3});
        Ars(Ars(:,4)<min(xx{4}),4) = min(xx{4});
        Ars(Ars(:,5)<min(xx{5}),5) = min(xx{5});
        Ars(Ars(:,6)<min(xx{6}),6) = min(xx{6});

        Ars(Ars(:,1)>max(xx{1}),1) = max(xx{1});
        Ars(Ars(:,2)>max(xx{2}),2) = max(xx{2});
        Ars(Ars(:,3)>max(xx{3}),3) = max(xx{3});
        Ars(Ars(:,4)>max(xx{4}),4) = max(xx{4});
        Ars(Ars(:,5)>max(xx{5}),5) = max(xx{5});
        Ars(Ars(:,6)>max(xx{6}),6) = max(xx{6});

        ind(:,1,1)  = 1 + floor(round((Ars(:,1)-xx{1}(1))/h(1)*1e8)*1e-8);
        ind(:,2,1)  = 1 + floor(round((Ars(:,2)-xx{2}(1))/h(2)*1e8)*1e-8);
        ind(:,3,1)  = 1 + floor(round((Ars(:,3)-xx{3}(1))/h(3)*1e8)*1e-8);
        ind(:,4,1)  = 1 + floor(round((Ars(:,4)-xx{4}(1))/h(4)*1e8)*1e-8);
        ind(:,5,1)  = 1 + floor(round((Ars(:,5)-xx{5}(1))/h(5)*1e8)*1e-8);
        ind(:,6,1)  = 1 + floor(round((Ars(:,6)-xx{6}(1))/h(6)*1e8)*1e-8);

        ind(:,1,2) = 1 +  ceil(round((Ars(:,1)-xx{1}(1))/h(1)*1e8)*1e-8);
        ind(:,2,2) = 1 +  ceil(round((Ars(:,2)-xx{2}(1))/h(2)*1e8)*1e-8);
        ind(:,3,2) = 1 +  ceil(round((Ars(:,3)-xx{3}(1))/h(3)*1e8)*1e-8);
        ind(:,4,2) = 1 +  ceil(round((Ars(:,4)-xx{4}(1))/h(4)*1e8)*1e-8);
        ind(:,5,2) = 1 +  ceil(round((Ars(:,5)-xx{5}(1))/h(5)*1e8)*1e-8);
        ind(:,6,2) = 1 +  ceil(round((Ars(:,6)-xx{6}(1))/h(6)*1e8)*1e-8);

        yy(:,1,1,1,1,1,1)    = YY(dpm_sub2ind(size(YY),ind(:,1,1),ind(:,2,1),ind(:,3,1),ind(:,4,1),ind(:,5,1),ind(:,6,1)));
        yy(:,2,1,1,1,1,1)    = YY(dpm_sub2ind(size(YY),ind(:,1,2),ind(:,2,1),ind(:,3,1),ind(:,4,1),ind(:,5,1),ind(:,6,1)));
        yy(:,1,2,1,1,1,1)    = YY(dpm_sub2ind(size(YY),ind(:,1,1),ind(:,2,2),ind(:,3,1),ind(:,4,1),ind(:,5,1),ind(:,6,1)));
        yy(:,2,2,1,1,1,1)    = YY(dpm_sub2ind(size(YY),ind(:,1,2),ind(:,2,2),ind(:,3,1),ind(:,4,1),ind(:,5,1),ind(:,6,1)));
        yy(:,1,1,2,1,1,1)    = YY(dpm_sub2ind(size(YY),ind(:,1,1),ind(:,2,1),ind(:,3,2),ind(:,4,1),ind(:,5,1),ind(:,6,1)));
        yy(:,2,1,2,1,1,1)    = YY(dpm_sub2ind(size(YY),ind(:,1,2),ind(:,2,1),ind(:,3,2),ind(:,4,1),ind(:,5,1),ind(:,6,1)));
        yy(:,1,2,2,1,1,1)    = YY(dpm_sub2ind(size(YY),ind(:,1,1),ind(:,2,2),ind(:,3,2),ind(:,4,1),ind(:,5,1),ind(:,6,1)));
        yy(:,2,2,2,1,1,1)    = YY(dpm_sub2ind(size(YY),ind(:,1,2),ind(:,2,2),ind(:,3,2),ind(:,4,1),ind(:,5,1),ind(:,6,1)));
        yy(:,1,1,1,2,1,1)    = YY(dpm_sub2ind(size(YY),ind(:,1,1),ind(:,2,1),ind(:,3,1),ind(:,4,2),ind(:,5,1),ind(:,6,1)));
        yy(:,2,1,1,2,1,1)    = YY(dpm_sub2ind(size(YY),ind(:,1,2),ind(:,2,1),ind(:,3,1),ind(:,4,2),ind(:,5,1),ind(:,6,1)));
        yy(:,1,2,1,2,1,1)    = YY(dpm_sub2ind(size(YY),ind(:,1,1),ind(:,2,2),ind(:,3,1),ind(:,4,2),ind(:,5,1),ind(:,6,1)));
        yy(:,2,2,1,2,1,1)    = YY(dpm_sub2ind(size(YY),ind(:,1,2),ind(:,2,2),ind(:,3,1),ind(:,4,2),ind(:,5,1),ind(:,6,1)));
        yy(:,1,1,2,2,1,1)    = YY(dpm_sub2ind(size(YY),ind(:,1,1),ind(:,2,1),ind(:,3,2),ind(:,4,2),ind(:,5,1),ind(:,6,1)));
        yy(:,2,1,2,2,1,1)    = YY(dpm_sub2ind(size(YY),ind(:,1,2),ind(:,2,1),ind(:,3,2),ind(:,4,2),ind(:,5,1),ind(:,6,1)));
        yy(:,1,2,2,2,1,1)    = YY(dpm_sub2ind(size(YY),ind(:,1,1),ind(:,2,2),ind(:,3,2),ind(:,4,2),ind(:,5,1),ind(:,6,1)));
        yy(:,2,2,2,2,1,1)    = YY(dpm_sub2ind(size(YY),ind(:,1,2),ind(:,2,2),ind(:,3,2),ind(:,4,2),ind(:,5,1),ind(:,6,1)));
        yy(:,1,1,1,1,2,1)    = YY(dpm_sub2ind(size(YY),ind(:,1,1),ind(:,2,1),ind(:,3,1),ind(:,4,1),ind(:,5,2),ind(:,6,1)));
        yy(:,2,1,1,1,2,1)    = YY(dpm_sub2ind(size(YY),ind(:,1,2),ind(:,2,1),ind(:,3,1),ind(:,4,1),ind(:,5,2),ind(:,6,1)));
        yy(:,1,2,1,1,2,1)    = YY(dpm_sub2ind(size(YY),ind(:,1,1),ind(:,2,2),ind(:,3,1),ind(:,4,1),ind(:,5,2),ind(:,6,1)));
        yy(:,2,2,1,1,2,1)    = YY(dpm_sub2ind(size(YY),ind(:,1,2),ind(:,2,2),ind(:,3,1),ind(:,4,1),ind(:,5,2),ind(:,6,1)));
        yy(:,1,1,2,1,2,1)    = YY(dpm_sub2ind(size(YY),ind(:,1,1),ind(:,2,1),ind(:,3,2),ind(:,4,1),ind(:,5,2),ind(:,6,1)));
        yy(:,2,1,2,1,2,1)    = YY(dpm_sub2ind(size(YY),ind(:,1,2),ind(:,2,1),ind(:,3,2),ind(:,4,1),ind(:,5,2),ind(:,6,1)));
        yy(:,1,2,2,1,2,1)    = YY(dpm_sub2ind(size(YY),ind(:,1,1),ind(:,2,2),ind(:,3,2),ind(:,4,1),ind(:,5,2),ind(:,6,1)));
        yy(:,2,2,2,1,2,1)    = YY(dpm_sub2ind(size(YY),ind(:,1,2),ind(:,2,2),ind(:,3,2),ind(:,4,1),ind(:,5,2),ind(:,6,1)));
        yy(:,1,1,1,2,2,1)    = YY(dpm_sub2ind(size(YY),ind(:,1,1),ind(:,2,1),ind(:,3,1),ind(:,4,2),ind(:,5,2),ind(:,6,1)));
        yy(:,2,1,1,2,2,1)    = YY(dpm_sub2ind(size(YY),ind(:,1,2),ind(:,2,1),ind(:,3,1),ind(:,4,2),ind(:,5,2),ind(:,6,1)));
        yy(:,1,2,1,2,2,1)    = YY(dpm_sub2ind(size(YY),ind(:,1,1),ind(:,2,2),ind(:,3,1),ind(:,4,2),ind(:,5,2),ind(:,6,1)));
        yy(:,2,2,1,2,2,1)    = YY(dpm_sub2ind(size(YY),ind(:,1,2),ind(:,2,2),ind(:,3,1),ind(:,4,2),ind(:,5,2),ind(:,6,1)));
        yy(:,1,1,2,2,2,1)    = YY(dpm_sub2ind(size(YY),ind(:,1,1),ind(:,2,1),ind(:,3,2),ind(:,4,2),ind(:,5,2),ind(:,6,1)));
        yy(:,2,1,2,2,2,1)    = YY(dpm_sub2ind(size(YY),ind(:,1,2),ind(:,2,1),ind(:,3,2),ind(:,4,2),ind(:,5,2),ind(:,6,1)));
        yy(:,1,2,2,2,2,1)    = YY(dpm_sub2ind(size(YY),ind(:,1,1),ind(:,2,2),ind(:,3,2),ind(:,4,2),ind(:,5,2),ind(:,6,1)));
        yy(:,2,2,2,2,2,1)    = YY(dpm_sub2ind(size(YY),ind(:,1,2),ind(:,2,2),ind(:,3,2),ind(:,4,2),ind(:,5,2),ind(:,6,1)));
        yy(:,1,1,1,1,1,2)    = YY(dpm_sub2ind(size(YY),ind(:,1,1),ind(:,2,1),ind(:,3,1),ind(:,4,1),ind(:,5,1),ind(:,6,2)));
        yy(:,2,1,1,1,1,2)    = YY(dpm_sub2ind(size(YY),ind(:,1,2),ind(:,2,1),ind(:,3,1),ind(:,4,1),ind(:,5,1),ind(:,6,2)));
        yy(:,1,2,1,1,1,2)    = YY(dpm_sub2ind(size(YY),ind(:,1,1),ind(:,2,2),ind(:,3,1),ind(:,4,1),ind(:,5,1),ind(:,6,2)));
        yy(:,2,2,1,1,1,2)    = YY(dpm_sub2ind(size(YY),ind(:,1,2),ind(:,2,2),ind(:,3,1),ind(:,4,1),ind(:,5,1),ind(:,6,2)));
        yy(:,1,1,2,1,1,2)    = YY(dpm_sub2ind(size(YY),ind(:,1,1),ind(:,2,1),ind(:,3,2),ind(:,4,1),ind(:,5,1),ind(:,6,2)));
        yy(:,2,1,2,1,1,2)    = YY(dpm_sub2ind(size(YY),ind(:,1,2),ind(:,2,1),ind(:,3,2),ind(:,4,1),ind(:,5,1),ind(:,6,2)));
        yy(:,1,2,2,1,1,2)    = YY(dpm_sub2ind(size(YY),ind(:,1,1),ind(:,2,2),ind(:,3,2),ind(:,4,1),ind(:,5,1),ind(:,6,2)));
        yy(:,2,2,2,1,1,2)    = YY(dpm_sub2ind(size(YY),ind(:,1,2),ind(:,2,2),ind(:,3,2),ind(:,4,1),ind(:,5,1),ind(:,6,2)));
        yy(:,1,1,1,2,1,2)    = YY(dpm_sub2ind(size(YY),ind(:,1,1),ind(:,2,1),ind(:,3,1),ind(:,4,2),ind(:,5,1),ind(:,6,2)));
        yy(:,2,1,1,2,1,2)    = YY(dpm_sub2ind(size(YY),ind(:,1,2),ind(:,2,1),ind(:,3,1),ind(:,4,2),ind(:,5,1),ind(:,6,2)));
        yy(:,1,2,1,2,1,2)    = YY(dpm_sub2ind(size(YY),ind(:,1,1),ind(:,2,2),ind(:,3,1),ind(:,4,2),ind(:,5,1),ind(:,6,2)));
        yy(:,2,2,1,2,1,2)    = YY(dpm_sub2ind(size(YY),ind(:,1,2),ind(:,2,2),ind(:,3,1),ind(:,4,2),ind(:,5,1),ind(:,6,2)));
        yy(:,1,1,2,2,1,2)    = YY(dpm_sub2ind(size(YY),ind(:,1,1),ind(:,2,1),ind(:,3,2),ind(:,4,2),ind(:,5,1),ind(:,6,2)));
        yy(:,2,1,2,2,1,2)    = YY(dpm_sub2ind(size(YY),ind(:,1,2),ind(:,2,1),ind(:,3,2),ind(:,4,2),ind(:,5,1),ind(:,6,2)));
        yy(:,1,2,2,2,1,2)    = YY(dpm_sub2ind(size(YY),ind(:,1,1),ind(:,2,2),ind(:,3,2),ind(:,4,2),ind(:,5,1),ind(:,6,2)));
        yy(:,2,2,2,2,1,2)    = YY(dpm_sub2ind(size(YY),ind(:,1,2),ind(:,2,2),ind(:,3,2),ind(:,4,2),ind(:,5,1),ind(:,6,2)));
        yy(:,1,1,1,1,2,2)    = YY(dpm_sub2ind(size(YY),ind(:,1,1),ind(:,2,1),ind(:,3,1),ind(:,4,1),ind(:,5,2),ind(:,6,2)));
        yy(:,2,1,1,1,2,2)    = YY(dpm_sub2ind(size(YY),ind(:,1,2),ind(:,2,1),ind(:,3,1),ind(:,4,1),ind(:,5,2),ind(:,6,2)));
        yy(:,1,2,1,1,2,2)    = YY(dpm_sub2ind(size(YY),ind(:,1,1),ind(:,2,2),ind(:,3,1),ind(:,4,1),ind(:,5,2),ind(:,6,2)));
        yy(:,2,2,1,1,2,2)    = YY(dpm_sub2ind(size(YY),ind(:,1,2),ind(:,2,2),ind(:,3,1),ind(:,4,1),ind(:,5,2),ind(:,6,2)));
        yy(:,1,1,2,1,2,2)    = YY(dpm_sub2ind(size(YY),ind(:,1,1),ind(:,2,1),ind(:,3,2),ind(:,4,1),ind(:,5,2),ind(:,6,2)));
        yy(:,2,1,2,1,2,2)    = YY(dpm_sub2ind(size(YY),ind(:,1,2),ind(:,2,1),ind(:,3,2),ind(:,4,1),ind(:,5,2),ind(:,6,2)));
        yy(:,1,2,2,1,2,2)    = YY(dpm_sub2ind(size(YY),ind(:,1,1),ind(:,2,2),ind(:,3,2),ind(:,4,1),ind(:,5,2),ind(:,6,2)));
        yy(:,2,2,2,1,2,2)    = YY(dpm_sub2ind(size(YY),ind(:,1,2),ind(:,2,2),ind(:,3,2),ind(:,4,1),ind(:,5,2),ind(:,6,2)));
        yy(:,1,1,1,2,2,2)    = YY(dpm_sub2ind(size(YY),ind(:,1,1),ind(:,2,1),ind(:,3,1),ind(:,4,2),ind(:,5,2),ind(:,6,2)));
        yy(:,2,1,1,2,2,2)    = YY(dpm_sub2ind(size(YY),ind(:,1,2),ind(:,2,1),ind(:,3,1),ind(:,4,2),ind(:,5,2),ind(:,6,2)));
        yy(:,1,2,1,2,2,2)    = YY(dpm_sub2ind(size(YY),ind(:,1,1),ind(:,2,2),ind(:,3,1),ind(:,4,2),ind(:,5,2),ind(:,6,2)));
        yy(:,2,2,1,2,2,2)    = YY(dpm_sub2ind(size(YY),ind(:,1,2),ind(:,2,2),ind(:,3,1),ind(:,4,2),ind(:,5,2),ind(:,6,2)));
        yy(:,1,1,2,2,2,2)    = YY(dpm_sub2ind(size(YY),ind(:,1,1),ind(:,2,1),ind(:,3,2),ind(:,4,2),ind(:,5,2),ind(:,6,2)));
        yy(:,2,1,2,2,2,2)    = YY(dpm_sub2ind(size(YY),ind(:,1,2),ind(:,2,1),ind(:,3,2),ind(:,4,2),ind(:,5,2),ind(:,6,2)));
        yy(:,1,2,2,2,2,2)    = YY(dpm_sub2ind(size(YY),ind(:,1,1),ind(:,2,2),ind(:,3,2),ind(:,4,2),ind(:,5,2),ind(:,6,2)));
        yy(:,2,2,2,2,2,2)    = YY(dpm_sub2ind(size(YY),ind(:,1,2),ind(:,2,2),ind(:,3,2),ind(:,4,2),ind(:,5,2),ind(:,6,2)));
        
        da(:,1) = (Ars(:,1)-xx{1}(ind(:,1,1)))/h(1);
        da(:,2) = (Ars(:,2)-xx{2}(ind(:,2,1)))/h(2);
        da(:,3) = (Ars(:,3)-xx{3}(ind(:,3,1)))/h(3);
        da(:,4) = (Ars(:,4)-xx{4}(ind(:,4,1)))/h(4);
        da(:,5) = (Ars(:,5)-xx{5}(ind(:,5,1)))/h(5);
        da(:,6) = (Ars(:,6)-xx{6}(ind(:,6,1)))/h(6);

        v5(:,1,1,1,1,1) = yy(:,1,1,1,1,1,1) + da(:,6).*(yy(:,1,1,1,1,1,2) - yy(:,1,1,1,1,1,1));
        v5(:,2,1,1,1,1) = yy(:,2,1,1,1,1,1) + da(:,6).*(yy(:,2,1,1,1,1,2) - yy(:,2,1,1,1,1,1));
        v5(:,1,2,1,1,1) = yy(:,1,2,1,1,1,1) + da(:,6).*(yy(:,1,2,1,1,1,2) - yy(:,1,2,1,1,1,1));
        v5(:,2,2,1,1,1) = yy(:,2,2,1,1,1,1) + da(:,6).*(yy(:,2,2,1,1,1,2) - yy(:,2,2,1,1,1,1));
        v5(:,1,1,2,1,1) = yy(:,1,1,2,1,1,1) + da(:,6).*(yy(:,1,1,2,1,1,2) - yy(:,1,1,2,1,1,1));
        v5(:,2,1,2,1,1) = yy(:,2,1,2,1,1,1) + da(:,6).*(yy(:,2,1,2,1,1,2) - yy(:,2,1,2,1,1,1));
        v5(:,1,2,2,1,1) = yy(:,1,2,2,1,1,1) + da(:,6).*(yy(:,1,2,2,1,1,2) - yy(:,1,2,2,1,1,1));
        v5(:,2,2,2,1,1) = yy(:,2,2,2,1,1,1) + da(:,6).*(yy(:,2,2,2,1,1,2) - yy(:,2,2,2,1,1,1));
        v5(:,1,1,1,2,1) = yy(:,1,1,1,2,1,1) + da(:,6).*(yy(:,1,1,1,2,1,2) - yy(:,1,1,1,2,1,1));
        v5(:,2,1,1,2,1) = yy(:,2,1,1,2,1,1) + da(:,6).*(yy(:,2,1,1,2,1,2) - yy(:,2,1,1,2,1,1));
        v5(:,1,2,1,2,1) = yy(:,1,2,1,2,1,1) + da(:,6).*(yy(:,1,2,1,2,1,2) - yy(:,1,2,1,2,1,1));
        v5(:,2,2,1,2,1) = yy(:,2,2,1,2,1,1) + da(:,6).*(yy(:,2,2,1,2,1,2) - yy(:,2,2,1,2,1,1));
        v5(:,1,1,2,2,1) = yy(:,1,1,2,2,1,1) + da(:,6).*(yy(:,1,1,2,2,1,2) - yy(:,1,1,2,2,1,1));
        v5(:,2,1,2,2,1) = yy(:,2,1,2,2,1,1) + da(:,6).*(yy(:,2,1,2,2,1,2) - yy(:,2,1,2,2,1,1));
        v5(:,1,2,2,2,1) = yy(:,1,2,2,2,1,1) + da(:,6).*(yy(:,1,2,2,2,1,2) - yy(:,1,2,2,2,1,1));
        v5(:,2,2,2,2,1) = yy(:,2,2,2,2,1,1) + da(:,6).*(yy(:,2,2,2,2,1,2) - yy(:,2,2,2,2,1,1));
        v5(:,1,1,1,1,2) = yy(:,1,1,1,1,2,1) + da(:,6).*(yy(:,1,1,1,1,2,2) - yy(:,1,1,1,1,2,1));
        v5(:,2,1,1,1,2) = yy(:,2,1,1,1,2,1) + da(:,6).*(yy(:,2,1,1,1,2,2) - yy(:,2,1,1,1,2,1));
        v5(:,1,2,1,1,2) = yy(:,1,2,1,1,2,1) + da(:,6).*(yy(:,1,2,1,1,2,2) - yy(:,1,2,1,1,2,1));
        v5(:,2,2,1,1,2) = yy(:,2,2,1,1,2,1) + da(:,6).*(yy(:,2,2,1,1,2,2) - yy(:,2,2,1,1,2,1));
        v5(:,1,1,2,1,2) = yy(:,1,1,2,1,2,1) + da(:,6).*(yy(:,1,1,2,1,2,2) - yy(:,1,1,2,1,2,1));
        v5(:,2,1,2,1,2) = yy(:,2,1,2,1,2,1) + da(:,6).*(yy(:,2,1,2,1,2,2) - yy(:,2,1,2,1,2,1));
        v5(:,1,2,2,1,2) = yy(:,1,2,2,1,2,1) + da(:,6).*(yy(:,1,2,2,1,2,2) - yy(:,1,2,2,1,2,1));
        v5(:,2,2,2,1,2) = yy(:,2,2,2,1,2,1) + da(:,6).*(yy(:,2,2,2,1,2,2) - yy(:,2,2,2,1,2,1));
        v5(:,1,1,1,2,2) = yy(:,1,1,1,2,2,1) + da(:,6).*(yy(:,1,1,1,2,2,2) - yy(:,1,1,1,2,2,1));
        v5(:,2,1,1,2,2) = yy(:,2,1,1,2,2,1) + da(:,6).*(yy(:,2,1,1,2,2,2) - yy(:,2,1,1,2,2,1));
        v5(:,1,2,1,2,2) = yy(:,1,2,1,2,2,1) + da(:,6).*(yy(:,1,2,1,2,2,2) - yy(:,1,2,1,2,2,1));
        v5(:,2,2,1,2,2) = yy(:,2,2,1,2,2,1) + da(:,6).*(yy(:,2,2,1,2,2,2) - yy(:,2,2,1,2,2,1));
        v5(:,1,1,2,2,2) = yy(:,1,1,2,2,2,1) + da(:,6).*(yy(:,1,1,2,2,2,2) - yy(:,1,1,2,2,2,1));
        v5(:,2,1,2,2,2) = yy(:,2,1,2,2,2,1) + da(:,6).*(yy(:,2,1,2,2,2,2) - yy(:,2,1,2,2,2,1));
        v5(:,1,2,2,2,2) = yy(:,1,2,2,2,2,1) + da(:,6).*(yy(:,1,2,2,2,2,2) - yy(:,1,2,2,2,2,1));
        v5(:,2,2,2,2,2) = yy(:,2,2,2,2,2,1) + da(:,6).*(yy(:,2,2,2,2,2,2) - yy(:,2,2,2,2,2,1));

        v4(:,1,1,1,1) = v5(:,1,1,1,1,1) + da(:,5).*(v5(:,1,1,1,1,2) - v5(:,1,1,1,1,1));
        v4(:,2,1,1,1) = v5(:,2,1,1,1,1) + da(:,5).*(v5(:,2,1,1,1,2) - v5(:,2,1,1,1,1));
        v4(:,1,2,1,1) = v5(:,1,2,1,1,1) + da(:,5).*(v5(:,1,2,1,1,2) - v5(:,1,2,1,1,1));
        v4(:,2,2,1,1) = v5(:,2,2,1,1,1) + da(:,5).*(v5(:,2,2,1,1,2) - v5(:,2,2,1,1,1));
        v4(:,1,1,2,1) = v5(:,1,1,2,1,1) + da(:,5).*(v5(:,1,1,2,1,2) - v5(:,1,1,2,1,1));
        v4(:,2,1,2,1) = v5(:,2,1,2,1,1) + da(:,5).*(v5(:,2,1,2,1,2) - v5(:,2,1,2,1,1));
        v4(:,1,2,2,1) = v5(:,1,2,2,1,1) + da(:,5).*(v5(:,1,2,2,1,2) - v5(:,1,2,2,1,1));
        v4(:,2,2,2,1) = v5(:,2,2,2,1,1) + da(:,5).*(v5(:,2,2,2,1,2) - v5(:,2,2,2,1,1));
        v4(:,1,1,1,2) = v5(:,1,1,1,2,1) + da(:,5).*(v5(:,1,1,1,2,2) - v5(:,1,1,1,2,1));
        v4(:,2,1,1,2) = v5(:,2,1,1,2,1) + da(:,5).*(v5(:,2,1,1,2,2) - v5(:,2,1,1,2,1));
        v4(:,1,2,1,2) = v5(:,1,2,1,2,1) + da(:,5).*(v5(:,1,2,1,2,2) - v5(:,1,2,1,2,1));
        v4(:,2,2,1,2) = v5(:,2,2,1,2,1) + da(:,5).*(v5(:,2,2,1,2,2) - v5(:,2,2,1,2,1));
        v4(:,1,1,2,2) = v5(:,1,1,2,2,1) + da(:,5).*(v5(:,1,1,2,2,2) - v5(:,1,1,2,2,1));
        v4(:,2,1,2,2) = v5(:,2,1,2,2,1) + da(:,5).*(v5(:,2,1,2,2,2) - v5(:,2,1,2,2,1));
        v4(:,1,2,2,2) = v5(:,1,2,2,2,1) + da(:,5).*(v5(:,1,2,2,2,2) - v5(:,1,2,2,2,1));
        v4(:,2,2,2,2) = v5(:,2,2,2,2,1) + da(:,5).*(v5(:,2,2,2,2,2) - v5(:,2,2,2,2,1));

        
        v3(:,1,1,1) = v4(:,1,1,1,1) + da(:,4).*(v4(:,1,1,1,2) - v4(:,1,1,1,1));
        v3(:,2,1,1) = v4(:,2,1,1,1) + da(:,4).*(v4(:,2,1,1,2) - v4(:,2,1,1,1));
        v3(:,1,2,1) = v4(:,1,2,1,1) + da(:,4).*(v4(:,1,2,1,2) - v4(:,1,2,1,1));
        v3(:,2,2,1) = v4(:,2,2,1,1) + da(:,4).*(v4(:,2,2,1,2) - v4(:,2,2,1,1));
        v3(:,1,1,2) = v4(:,1,1,2,1) + da(:,4).*(v4(:,1,1,2,2) - v4(:,1,1,2,1));
        v3(:,2,1,2) = v4(:,2,1,2,1) + da(:,4).*(v4(:,2,1,2,2) - v4(:,2,1,2,1));
        v3(:,1,2,2) = v4(:,1,2,2,1) + da(:,4).*(v4(:,1,2,2,2) - v4(:,1,2,2,1));
        v3(:,2,2,2) = v4(:,2,2,2,1) + da(:,4).*(v4(:,2,2,2,2) - v4(:,2,2,2,1));
        
        v2(:,1,1) = v3(:,1,1,1) + da(:,3).*(v3(:,1,1,2) - v3(:,1,1,1));
        v2(:,2,1) = v3(:,2,1,1) + da(:,3).*(v3(:,2,1,2) - v3(:,2,1,1));
        v2(:,1,2) = v3(:,1,2,1) + da(:,3).*(v3(:,1,2,2) - v3(:,1,2,1));
        v2(:,2,2) = v3(:,2,2,1) + da(:,3).*(v3(:,2,2,2) - v3(:,2,2,1)); 
        
        v1(:,1) = v2(:,1,1) + da(:,2).*(v2(:,1,2) - v2(:,1,1));
        v1(:,2) = v2(:,2,1) + da(:,2).*(v2(:,2,2) - v2(:,2,1));

        yi = da(:,1).*(v1(:,2) - v1(:,1)) + v1(:,1);            
        
        y = reshape(yi,size(A1));
    otherwise
        error('DPM:Internal','Too many states or inputs: contact the author of DPM')
end
