CREATE TABLE tmp_nopp
SELECT target,method,gdt 
FROM benchmark 
WHERE cct=4 and max=10 and phipsi='-' and omega='-' and gdt!=0 
ORDER BY target,method;

CREATE TABLE tmp_pp20
SELECT target,method,gdt 
FROM benchmark 
WHERE cct=4 and max=10 and phipsi='pp20' and omega='-' and gdt!=0 
ORDER BY target,method;

CREATE TABLE tmp_pp20mm0
SELECT target,method,gdt 
FROM benchmark 
WHERE cct=4 and max=10 and phipsi='pp20mm0' and omega='-' and gdt!=0 
ORDER BY target,method;

CREATE TABLE tmp_pp30 
SELECT target,method,gdt 
FROM benchmark 
WHERE cct=4 and max=10 and phipsi='pp30' and omega='-' and gdt!=0 
ORDER BY target,method;

CREATE TABLE tmp_pp20mm2_om178-5 
SELECT target,method,gdt 
FROM benchmark 
WHERE cct=4 and max=10 and phipsi='pp20mm2' and omega='om178-5' and gdt!=0 
ORDER BY target,method;

CREATE TABLE tmp_targets_sorted 
SELECT substring(target,1,5) as target, avg(gdt_ts) as gdt 
FROM results_raw 
GROUP by substring(target,1,5) 
ORDER BY gdt desc;


-- to plot nopp vs pp20 vs pp30 vs pp20mm0
CREATE TABLE tmp 
SELECT a.target as target,a.method as method, a.gdt as gdt_nopp, b.gdt as gdt_pp20, c.gdt as gdt_pp30, d.gdt as gdt_pp20mm0
FROM tmp_nopp as a 
INNER JOIN tmp_pp20 as b 
INNER JOIN tmp_pp30 as c 
INNER JOIN tmp_pp20mm0 as d
INNER JOIN tmp_targets_sorted as e
ON (a.target=b.target AND a.target=c.target AND a.target=d.target AND a.method=b.method AND a.method=c.method AND a.method=d.method AND a.target=e.target) 
ORDER BY e.gdt desc;

-- same as above only for blast 
CREATE TABLE tmp_bla
SELECT a.target as target,a.method as method, a.gdt as gdt_nopp, b.gdt as gdt_pp20, c.gdt as gdt_pp30, d.gdt as gdt_pp20mm0
FROM tmp_nopp as a 
INNER JOIN tmp_pp20 as b 
INNER JOIN tmp_pp30 as c
INNER JOIN tmp_pp20mm0 as d 
INNER JOIN tmp_targets_sorted as e
ON (a.target=b.target AND a.target=c.target AND a.target=d.target AND a.method=b.method AND a.method=c.method AND a.method=d.method AND a.target=e.target) 
WHERE a.method='bla'
ORDER BY e.gdt desc;

-- same only for psi-blast
CREATE TABLE tmp_pb3
SELECT a.target as target,a.method as method, a.gdt as gdt_nopp, b.gdt as gdt_pp20, c.gdt as gdt_pp30, d.gdt as gdt_pp20mm0
FROM tmp_nopp as a 
INNER JOIN tmp_pp20 as b 
INNER JOIN tmp_pp30 as c 
INNER JOIN tmp_pp20mm0 as d
INNER JOIN tmp_targets_sorted as e
ON (a.target=b.target AND a.target=c.target AND a.target=d.target AND a.method=b.method AND a.method=c.method AND a.method=d.method AND a.target=e.target) 
WHERE a.method='pb3'
ORDER BY e.gdt desc;

-- same only for gtg
CREATE TABLE tmp_gtg
SELECT a.target as target,a.method as method, a.gdt as gdt_nopp, b.gdt as gdt_pp20, c.gdt as gdt_pp30, d.gdt as gdt_pp20mm0
FROM tmp_nopp as a 
INNER JOIN tmp_pp20 as b 
INNER JOIN tmp_pp30 as c 
INNER JOIN tmp_pp20mm0 as d
INNER JOIN tmp_targets_sorted as e
ON (a.target=b.target AND a.target=c.target AND a.target=d.target AND a.method=b.method AND a.method=c.method AND a.method=d.method AND a.target=e.target) 
WHERE a.method='gtg'
ORDER BY e.gdt desc;

-- blast vs psi-blast vs gtg (using pp20, max10, cct4)
SELECT a.target as target, a.method as method, a.gdt_pp20 as gdt_bla, b.gdt_pp20 as gdt_pb3, c.gdt_pp20 as gdt_gtg
FROM tmp AS a
JOIN tmp AS b
JOIN tmp AS c
ON (a.target=b.target AND a.target=c.target)
WHERE a.method='bla' AND b.method='pb3' AND c.method='gtg';

-- psi-blast vs gtg (using pp20, max10, cct4)
SELECT b.target as target, b.method as method, b.gdt_pp20 as gdt_pb3, c.gdt_pp20 as gdt_gtg
FROM tmp AS b
JOIN tmp AS c
ON (b.target=c.target)
WHERE b.method='pb3' AND c.method='gtg';


-- nopp vs pp20mm0 vs pp20mm2_om178-5 (blast)
SELECT a.target as target,a.method as method, a.gdt as gdt_nopp, b.gdt as gdt_pp20, c.gdt as gdt_pp20mm0, d.gdt as gdt_pp20mm2_om178_5
FROM tmp_nopp as a 
INNER JOIN tmp_pp20 as b  
INNER JOIN tmp_pp20mm0 as c
INNER JOIN tmp_pp20mm2_om178_5 as d
INNER JOIN tmp_targets_sorted as e
ON (a.target=b.target AND a.target=c.target AND a.target=d.target AND a.method=b.method AND a.method=c.method AND a.method=d.method AND a.target=e.target) 
WHERE a.method='bla'
ORDER BY e.gdt desc;

-- nopp vs pp20mm0 vs pp20mm2_om178-5 (psi-blast)
SELECT a.target as target,a.method as method, a.gdt as gdt_nopp, b.gdt as gdt_pp20, c.gdt as gdt_pp20mm0, d.gdt as gdt_pp20mm2_om178_5
FROM tmp_nopp as a 
INNER JOIN tmp_pp20 as b  
INNER JOIN tmp_pp20mm0 as c
INNER JOIN tmp_pp20mm2_om178_5 as d
INNER JOIN tmp_targets_sorted as e
ON (a.target=b.target AND a.target=c.target AND a.target=d.target AND a.method=b.method AND a.method=c.method AND a.method=d.method AND a.target=e.target) 
WHERE a.method='pb3'
ORDER BY e.gdt desc;

-- nopp vs pp20mm0 vs pp20mm2_om178-5 (gtg)
SELECT a.target as target,a.method as method, a.gdt as gdt_nopp, b.gdt as gdt_pp20, c.gdt as gdt_pp20mm0, d.gdt as gdt_pp20mm2_om178_5
FROM tmp_nopp as a 
INNER JOIN tmp_pp20 as b  
INNER JOIN tmp_pp20mm0 as c
INNER JOIN tmp_pp20mm2_om178_5 as d
INNER JOIN tmp_targets_sorted as e
ON (a.target=b.target AND a.target=c.target AND a.target=d.target AND a.method=b.method AND a.method=c.method AND a.method=d.method AND a.target=e.target) 
WHERE a.method='gtg'
ORDER BY e.gdt desc;
