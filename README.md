# OsciloscopeGenerateRAWEMZ

Octave file for generate .RAWEMZ from .xlsx files

- New Folder 1 - 3 คือไฟล์สำหรับลองสร้าง rawemz ฟังก์ชัน

- New Folder 4 คือฟังก์ชันสำหรับ generate rawemz 64 channel จากไฟล์ xlsx 64 ไฟล์
  ถ้าใส่ไม่ถึง64ไฟล์ มันจะ error

- New Folder 5 คือฟังก์ชันสำหรับ generate rawemz 64 channel เหมือนกับ New Folder 4
  แต่ฟังก์ชันนี้จะเช็คว่าไฟล์ที่รับมา ถ้าไม่ถึง 64 ไฟล์ก็จะเติม channel ที่ขาดให้ (เติม0)
