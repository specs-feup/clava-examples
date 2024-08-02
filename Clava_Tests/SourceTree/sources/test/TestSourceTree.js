import lara.Check;

import clava.Clava;

/**
 * @test
 */
aspectdef TestsSourceTree
	var files = [];

	for(var $file of Clava.getProgram().children) {
		var filepath = $file.relativeFilepath;
		if($file.sourceFoldername !== undefined) {
			filepath = $file.sourceFoldername + "/" + filepath;
		}

		files.push(filepath);
	}

	var result = files.sort().join();
	var expected = "nested/srcFolderExtra2/lib/lib4.c,nested/srcFolderExtra2/lib/lib_nested/lib5.c,nested/srcFolderExtra2/test8.c,src/lib/lib1.c,src/test1.c,src/test2.c,src2/lib/lib2.c,src2/lib/lib2/lib6.c,src2/test5.c,src2/test6.c,srcFolderExtra/lib/lib3.c,srcFolderExtra/test7.c,test3.c,test4.c,test9.c";
	
	Check.strings(result, expected);
	
end




