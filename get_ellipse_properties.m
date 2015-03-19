function properties = get_ellipse_properties(bords)

regions = bwlabel(logical(1-bords), 4);
props = regionprops(regions,'MajorAxisLength','MinorAxisLength','Orientation');
major = [props.MajorAxisLength];
minor = [props.MinorAxisLength];
orient = [props.Orientation];

%gets rid of the background
major = removerows(major(:), 1);
minor = removerows(minor(:), 1);
orient = removerows(orient(:), 1);

properties = [major(:) minor(:) orient(:)];